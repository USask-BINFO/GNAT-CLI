#!/usr/bin/env python3

"""
Class to faciliate label file manipulation.
"""

from decimal import Decimal
import colorsys
import random

class LabelObject:
    def __init__(self):
        self.accessions = {}
        self.is_new = False
        self.num_classes = 0
        self.classes = None
        self.header = None

    def from_file(self, label_file):
        """
        Parse the label file and store the labels in the object.
        """
        try:
            with open(label_file, 'r') as f:
                label_lines = f.readlines()
        except FileNotFoundError:
            raise FileNotFoundError(f"Label file {label_file} not found. Please check the path and try again.")
        
        if len(label_lines) > 0:
            self.is_new = self.is_new_format(label_lines)
            if self.is_new:
                self.__from_new_format(label_lines)
            else:
                self.__from_old_format(label_lines)

    def is_new_format(self, label_lines):
        """
        Check if the label file is in the new format or the old format.
        """
        if label_lines[0].startswith("#"):
            return True
        else:
            if len(label_lines[0].strip().split()) == 3:
                return False
            else:
                raise ValueError("Label file is not in the new or old format. Please check the file and try again. (First label line does not start with #, but label lines are also not 3 tab-delimited items.)")

    def __from_new_format(self, label_lines):
        """
        Parse the label file in the new format and build the label object.

        New format starts with a header line beginning with #, is comma-delimited, and has the format:
            accession,label,e_value,group1_name,group1_colours,group2_name,group2_colours,...,groupN_name,groupN_colours
        Where all items after [2] are optional. It is also possible for some groups to be empty.
        """
        self.header = label_lines[0].strip()
        split_header = self.header.split(',')
        self.num_classes = (len(split_header) - 3) // 2
        self.classes = [{'label_name': split_header[i], 'label_map': {}} for i in range(3, len(split_header), 2)]

        for line in label_lines[1:]:
            split_line = line.strip().split(',')
            self.add_line(split_line)

    def __from_old_format(self, label_lines):
        """
        Parse the label file in the old format and build the label object.

        Old format is tab-delimited and has the format:
            accession   label   e_value
        """
        self.header = None
        self.num_classes = 0
        self.classes = None

        for line in label_lines:
            split_line = line.strip().split()
            self.add_line(split_line)

    def add_line(self, split_line):
        """
        Add a line from the label file to the label object.
        """
        accession, label, e_value = split_line[:3]
        replace_dupe = False

        # print(f"Adding line: {split_line}")

        if self.is_new:
            # Add line in new format
            if accession in self.accessions:
                if len(split_line[3:]) != 2 * self.num_classes:
                    # Pad with blanks
                    split_line += ["" for _ in range(2 * self.num_classes - len(split_line[3:]))]

                # Accession already exists, update label and e_value with the lower of the two e_values
                if Decimal(e_value) < Decimal(self.accessions[accession][1]):
                    replace_dupe = True
                    old_classes = self.accessions[accession][2]
                    self.accessions[accession] = [label, e_value, split_line[3:]]

                else:
                    # New dupilcate has a larger e_value so we can just return
                    return
            else:
                # New accession
                self.accessions[accession] = [label, e_value, split_line[3:]]

            if replace_dupe:
                # TODO: SHOULD I DO THIS? SHOULD I DELETE THE OLD CLASS LABELS??? THERE ISN'T A NEW CLASS LABEL... nevermind looks like I do add it later
                # Remove accession from old class label
                for j in range(0, len(old_classes), 2):
                    class_idx = j // 2
                    class_label = old_classes[j]

                    while accession in self.classes[class_idx]["label_map"][class_label]["accessions"]:
                        self.classes[class_idx]["label_map"][class_label]["accessions"].remove(accession)

            for j in range(3, len(split_line), 2):
                class_idx = (j - 3) // 2

                class_label = split_line[j]
                class_colour = split_line[j+1]

                if class_label not in self.classes[class_idx]["label_map"]:
                    self.classes[class_idx]["label_map"][class_label] = {"accessions": [accession], "colour": class_colour}

                else:
                    # Class label already exists, add accession to it (and make sure the colour is the same)
                    if class_colour != self.classes[class_idx]["label_map"][class_label]["colour"]:
                        print(f"Warning: Class label '{class_label}' has multiple colours in label file (new colour: '{class_colour}'). Using the first one in the label file: {self.classes[class_idx]['label_map'][class_label]['colour']}.")

                    # Double check that the accession isn't already in the class label (if it is, something is wrong)
                    if accession in self.classes[class_idx]["label_map"][class_label]["accessions"]:
                        raise ValueError(f"Accession {accession} is already in class label '{class_label}' at the second check. Something is very wrong.")
                    
                    self.classes[class_idx]["label_map"][class_label]["accessions"].append(accession)
        
        else:
            # Add line in old format
            if accession in self.accessions:
                # Accession already exists, update label and e_value with the lower of the two e_values
                if Decimal(e_value) < Decimal(self.accessions[accession][1]):
                    self.accessions[accession] = (label, e_value)
                else:
                    return
            else:
                self.accessions[accession] = (label, e_value)
            
    def add_parsed_line(self, split_line):
        """
        Add a line from a parsed PSI-BLAST output file to the label object.
        :param split_line: The split line from the PSI-BLAST output file (lines are originially tab-delimited).
            Format: [0]gbff_name, [1]match_accession, [2]e_value, [3]query_accession, [4]pid, [5]query_len, [6]match_len, [7]match_len/query_len
        """
        accession = split_line[1]
        label = split_line[3]
        e_value = split_line[2]

        self.add_line([accession, label, e_value])

    def merge_parsed_psi_file(self, parsed_location):
        """
        Merge parsed PSI-BLAST output lines with the label object.
        :param parsed_location: The location of the parsed PSI-BLAST output file.
        """
        try:
            with open(parsed_location, 'r') as f:
                parsed_lines = f.readlines()
        except FileNotFoundError:
            raise FileNotFoundError(f"Parsed PSI-BLAST output file {parsed_location} not found. Please check the path and try again.")
        
        for line in parsed_lines:
            self.add_parsed_line(line.split())

    def __str__(self, fmt="new"):
        """
        Convert the label object to a string of label file lines.
        """
        fmt = fmt == "new"
        if fmt != self.is_new:
            self.switch_format()
            
        if self.is_new:
            # New format
            label_str = f"{self.header}\n"
            for accession in self.accessions:
                label_str += f"{accession},{self.accessions[accession][0]},{self.accessions[accession][1]}"
                for class_idx in range(len(self.classes)):
                    add_flag = False
                    for class_label in self.classes[class_idx]["label_map"]:
                        if accession in self.classes[class_idx]["label_map"][class_label]["accessions"]:
                            label_str += f",{class_label},{self.classes[class_idx]['label_map'][class_label]['colour']}"
                            add_flag = True
                            break
                    if not add_flag:
                        label_str += ",,"
                label_str += "\n"
        else:
            # Old format
            label_str = ""
            for accession in self.accessions:
                label_str += f"{accession}\t{self.accessions[accession][0]}\t{self.accessions[accession][1]}\n"

        return label_str
    
    def propagate_class_labels(self, text_label, class_level, class_label, class_colour):
        """
        Assign a given class label/colour to all accessions with a given old (text) label.
        :param text_label: The text label to match. (E.g., "AcrIF7_protein")
        :param class_level: The index of the class level to assign the class label to. (E.g., 0)
        :param class_label: The class label to assign. (E.g., "anti-CRISPR")
        :param class_colour: The hex colour to assign to the class label. (E.g., "#FF0000")
        """
        if class_level >= len(self.classes):
            raise ValueError(f"Class level index {class_level} is out of range. There are only {len(self.classes)} class levels.")

        if self.is_new:
            # New format
            for accession in self.accessions:
                if self.accessions[accession][0] == text_label:
                    if class_label not in self.classes[class_level]["label_map"]:
                        self.classes[class_level]["label_map"][class_label] = {"accessions": [accession], "colour": class_colour}
                    else:
                        self.classes[class_level]["label_map"][class_label]["accessions"].append(accession)
        else:
            # Old format
            print(f"Warning: Propagating class labels is not supported for old format label files. Skipping propagation for text label '{text_label}'.")

    def add_class_level(self, class_label):
        """
        Add a new class level to the label object.
        :param class_label: The name of the new class level.
        """
        if self.is_new:
            # New format
            self.num_classes += 1
            self.classes.append({"label_name": class_label, "label_map": {}})
            self.header += f",{class_label},{class_label}_colour"

            for accession in self.accessions:
                self.accessions[accession][2].append(",,")

        else:
            # Old format
            print(f"Warning: Adding class levels is not supported for old format label files. Skipping addition of class level '{class_label}'.")

    def check_class_colours(self):
        """
        Make sure all class labels have consistent colouring/a colour is assigned rather than 'None'.
        """
        existing_colour_map = {}
        if self.is_new:
            # New format
            for class_idx in range(len(self.classes)):
                for class_label in self.classes[class_idx]["label_map"]:
                    if self.classes[class_idx]["label_map"][class_label]["colour"] is None:
                        print(f"Warning: Class label '{class_label}' has no colour assigned. Assigning a random colour.")
                        self.classes[class_idx]["label_map"][class_label]["colour"] = "#000000"
                    else:
                        colour = self.classes[class_idx]["label_map"][class_label]["colour"]
                        if colour in existing_colour_map:
                            if class_label not in existing_colour_map[colour]:
                                print(f"Warning: Colour '{colour}' is assigned to multiple class labels: '{existing_colour_map[colour]}' and '{class_label}'.")
                                existing_colour_map[colour].append(class_label)
                        else:
                            existing_colour_map[colour] = [class_label]
        else:
            # Old format
            print(f"Warning: Checking class colours is not supported for old format label files. Skipping colour check.")

    def colour_picker(self, n):
        """
        Generate n hex colours evenly spaced around the colour wheel.
        :param n: The number of colours to generate. (if n is 1, generates a random colour)
        """
        if n == 1:
            # Generate a random colour if n is 1
            hue = random.random() # random value 0..1
            r,g,b = colorsys.hsv_to_rgb(hue,1,1)
            return [f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"]
        
        colours = []
        for i in range(n):
            hue = i / n

            r,g,b = colorsys.hsv_to_rgb(hue,1,1)

            hex_colour = f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"

            colours.append(hex_colour)
        
        return colours

    def re_assign_colours(self, class_idx=None):
        """
        Re-assign colours to class labels.
        :param class_idx: The index of the class level to re-assign colours for. If None, re-assign colours for all class levels.
        """
        if self.is_new:
            # New format, re-assign the colours
            if class_idx is None:
                for i in range(len(self.classes)):
                    self.re_assign_colours(i)
            else:
                # Re-assign colours for a specific class level
                class_labels = list(self.classes[class_idx]["label_map"].keys())
                # removed_blanks = []
                for blank in ["", "None"]:
                    if blank in class_labels:
                        class_labels.remove(blank)
                        # removed_blanks.append(blank)
                colours = self.colour_picker(len(class_labels))

                for i in range(len(class_labels)):
                    self.classes[class_idx]["label_map"][class_labels[i]]["colour"] = colours[i]
        else:
            # Old format
            print(f"Warning: Re-assigning colours is not supported for old format label files. Skipping colour re-assignment.")

    def rename_label(self, old_label, new_label):
        if self.is_new:
            # New format, need to change the label in the accessions and classes
            for accession in self.accessions:
                if self.accessions[accession][0] == old_label:
                    self.accessions[accession][0] = new_label

            for class_idx in range(len(self.classes)):
                if old_label in self.classes[class_idx]["label_map"]:
                    self.classes[class_idx]["label_map"][new_label] = self.classes[class_idx]["label_map"].pop(old_label)
        else:
            # Old format
            for accession in self.accessions:
                if self.accessions[accession][0] == old_label:
                    self.accessions[accession][0] = new_label

    def remove_label(self, label):
        if self.is_new:
            # New format, need to remove the label from the accessions and classes
            for accession in self.accessions:
                if self.accessions[accession][0] == label:
                    del self.accessions[accession]

            for class_idx in range(len(self.classes)):
                if label in self.classes[class_idx]["label_map"]:
                    del self.classes[class_idx]["label_map"][label]
        else:
            # Old format
            for accession in self.accessions:
                if self.accessions[accession][0] == label:
                    del self.accessions[accession]

    def relabel_accession(self, accession, new_label, relabel_if_last=True):
        """
        :param relabel_if_last: If True, will relabel (new version) class labels with the new label if accession was the last one with the old label. If the new label already has a class label, the old one will be preserved. Default is False.
        """
        if self.is_new:
            # New format, need to change the label in the accessions
            if accession in self.accessions:
                old_label = self.accessions[accession][0]
                self.accessions[accession][0] = new_label

            # Check if that was the last accession with the old label
            old_label_count = 0
            for acc in self.accessions:
                if self.accessions[acc][0] == old_label:
                    old_label_count += 1

            if old_label_count == 0 and relabel_if_last:
                # No more accessions with the old label, relabel all class labels with the old label to the new label
                for class_idx in range(len(self.classes)):
                    if old_label in self.classes[class_idx]["label_map"]:
                        if new_label in self.classes[class_idx]["label_map"]:
                            print(f"Warning: Class label '{new_label}' already exists in class level {class_idx}. Skipping class label modification.")
                        else:
                            self.classes[class_idx]["label_map"][new_label] = self.classes[class_idx]["label_map"].pop(old_label)

        else:
            # Old format
            if accession in self.accessions:
                self.accessions[accession][0] = new_label

    def merge_label_file(self, new_label_file):
        """
        Merge a new label file with the current label object. Merged file will be in same format as old label object.
        :param new_label_file: The path to the new label file to merge.
        """
        # TODO DOUBLE CHECK THIS WORKS
        new_label_obj = LabelObject()
        new_label_obj.from_file(new_label_file)

        new_label_fmt = new_label_obj.is_new
        old_label_fmt = self.is_new

        if new_label_fmt != old_label_fmt:
            new_label_obj.switch_format()

        if old_label_fmt:
            # New format
            # Find how many new classes there are
            new_classes = [class_label["label_name"] for class_label in new_label_obj.classes]
            old_classes = [class_label["label_name"] for class_label in self.classes]

            index_map = {}
            counter = len(old_classes)
            for i, c in enumerate(new_classes):
                if c not in old_classes:
                    index_map[i] = counter
                    counter += 1
                else:
                    index_map[i] = old_classes.index(c)

            unique_new_classes = [(i, c) for i, c in enumerate(new_classes) if c not in old_classes]
            for i, c in unique_new_classes:
                # Add the new class level for each new class
                self.add_class_level(c)

            # Loop through the accessions to add them to the old label object
            for accession in new_label_obj.accessions:
                # Rebuild the line to add to the old label object
                new_label_line = [accession, new_label_obj.accessions[accession][0], new_label_obj.accessions[accession][1]]

                # Pad the line with blanks for the new total number of classes
                new_label_line += ["" for _ in range(2 * len(self.classes) - len(new_label_line[3:]))]

                line_new_classes = new_label_obj.accessions[accession][2]
                # Remap these new classes to the mapping
                for i in range(0, len(line_new_classes), 2):
                    print(line_new_classes)
                    j = i // 2
                    class_idx = index_map[j]
                    line_idx = 3 + 2 * class_idx
                    new_label_line[line_idx] = line_new_classes[i]
                    new_label_line[line_idx + 1] = self.colour_picker(1)[0]

                self.add_line(new_label_line)
            
            # Re-assign colours to the class levels
            self.re_assign_colours()

        else:
            # Old format
            for accession in new_label_obj.accessions:
                self.add_line([accession, new_label_obj.accessions[accession][0], new_label_obj.accessions[accession][1]])



    def switch_format(self):
        """
        Converts the label object to the opposite format. (Note that this will lose information if converting from new to old).
        """
        self.is_new = not self.is_new
        if self.is_new:
            # Convert to new format
            self.header = "#accession,label,e_value"
            self.num_classes = 0
            self.classes = []
            
            for accession in self.accessions:
                old_tuple = self.accessions[accession]
                self.accessions[accession] = [old_tuple[0], old_tuple[1], []]

        else:
            # Convert to old format
            self.header = None
            self.num_classes = 0
            self.classes = None

            for accession in self.accessions:
                old_list = self.accessions[accession]
                self.accessions[accession] = (old_list[0], old_list[1])



if __name__ == "__main__":
    # Testing
    print("Reading LF 1")
    label_obj_1 = LabelObject()
    label_obj_1.from_file("/home/njt694/gnat/utils/AcrVis/scripts/lf_testing/l1_new.label")
    print(label_obj_1)
    print(label_obj_1.classes)
    print(label_obj_1.accessions)

    print("\n\n")
    print("Reading LF 2")

    label_obj_2 = LabelObject()
    label_obj_2.from_file("/home/njt694/gnat/utils/AcrVis/scripts/lf_testing/l2_new.label")
    print(label_obj_2)
    print(label_obj_2.classes)
    print(label_obj_2.accessions)

    print("\n\n")
    print("Merging LF 2 into LF 1")

    label_obj_1.merge_label_file("/home/njt694/gnat/utils/AcrVis/scripts/lf_testing/l2_new.label")
    print(label_obj_1)
    print(label_obj_1.classes)
    print(label_obj_1.accessions)

    print("\n\n")
    print("Switching LF 1 format")

    label_obj_1.switch_format()
    print(label_obj_1)
    print(label_obj_1.classes)
    print(label_obj_1.accessions)

    print("\n\n")
    print("Reading LF 3 (old fmt)")

    label_obj_3 = LabelObject()
    label_obj_3.from_file("/home/njt694/gnat/utils/AcrVis/scripts/lf_testing/l3_old.label")
    print(label_obj_3)
    print(label_obj_3.classes)
    print(label_obj_3.accessions)

    print("\n\n")
    print("Switching LF 3 format")

    label_obj_3.switch_format()
    print(label_obj_3)
    print(label_obj_3.classes)
    print(label_obj_3.accessions)

    print("\n\n")
    print("Adding class level to LF 3")

    label_obj_3.add_class_level("new_class")
    print(label_obj_3)
    print(label_obj_3.classes)
    print(label_obj_3.accessions)

    print("\n\n")
    print("Adding class label to LF 3 Cas1f")

    label_obj_3.propagate_class_labels("Cas1f", 0, "new_class", "#FF0000")

    print(label_obj_3)
    print(label_obj_3.classes)
    print(label_obj_3.accessions)
    pass