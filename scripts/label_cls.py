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
            if len(label_lines[0].strip().split('\t')) == 3:
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
            split_line = line.strip().split('\t')
            self.add_line(split_line)

    def add_line(self, split_line):
        """
        Add a line from the label file to the label object.
        """
        accession, label, e_value = split_line[:3]
        replace_dupe = False

        if self.is_new:
            # Add line in new format
            if accession in self.accessions:
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

            for j in range(3, len(split_line, 2)):
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



    def __str__(self):
        """
        Convert the label object to a string of label file lines.
        """
        if self.is_new:
            # New format
            label_str = f"{self.header}\n"
            for accession in self.accessions:
                label_str += f"{accession},{self.accessions[accession][0]},{self.accessions[accession][1]}"
                for class_idx in range(len(self.classes)):
                    for class_label in self.classes[class_idx]["label_map"]:
                        if accession in self.classes[class_idx]["label_map"][class_label]["accessions"]:
                            label_str += f",{class_label},{self.classes[class_idx]['label_map'][class_label]['colour']}"
                        else:
                            label_str += ",,"
                label_str += "\n"
        else:
            # Old format
            label_str = ""
            for accession in self.accessions:
                label_str += f"{accession}\t{self.accessions[accession][0]}\t{self.accessions[accession][1]}\n"

        return label_str.rstrip() # Remove trailing newline
    
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
                    self.classes[i]["label_map"] = self.re_assign_colours(i)
            else:
                # Re-assign colours for a specific class level
                class_labels = list(self.classes[class_idx]["label_map"].keys())
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

if __name__ == "__main__":
    # Testing
    label_obj = LabelObject()

    print(label_obj.accessions)
    print(label_obj.classes)
    print(label_obj.header)
    print(label_obj.num_classes)
    pass