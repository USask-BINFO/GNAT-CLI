from os import path
from sortedcontainers import SortedList
from LabelRecord import LabelRecord

class LabelFile:
    """
    A wrapper to interact with a SortedList storing LabelRecord objects.
    This class provides methods to add, remove, and find records in the list.
    It also provides methods to bulk extract records and to merge two label files.

    The output file will be sorted by the accessions (first column), even if the 
    input file is not sorted.
    """

    def __init__(self, filename=None, max_file_size=2 * 1024 ** 3, ):
        """
        Initialize the LabelFile object.

        :param filename: The name of the file to read from or write to.
        :param max_file_size: The maximum size of the file in bytes. (default=2 GB)
        """
        self.filename = filename

        self.max_file_size = max_file_size

        # Even though I implemented __gt__ etc. in SortedRecord, for some reason explicitly defining the key is faster
        self.sorted_list = SortedList(key=lambda x: x.key) 
        
        self.header = None

        self.delim = ','

        if self.filename is not None:
            self.load_from_file(self.filename)

    def __str__(self):
        mystr = f"LabelFile({self.filename})\n"
        for line in self.sorted_list:
            mystr += str(line) + "\n"

        return mystr.rstrip()
    
    def __iter__(self):
        """
        Iterate over the sorted list.
        :return: An iterator over the sorted list.
        """
        return iter(self.sorted_list)

    def __human_readable_size__(self, size):
        """
        Convert a size in bytes to a human-readable format.
        :param size: Size in bytes.
        :return: Human-readable size string.
        """
        for unit in ['B', 'KB', 'MB', 'GB']:
            if size < 1024:
                return f"{size:.3f} {unit}"
            size /= 1024
        return f"{size:.3f} TB"
    
    def __verify_line_vs_header__(self, line):
        """
        Verify that the line matches the header format.
        :param line: The line to verify.
        :return: True if the line matches the header format, False otherwise.
        """
        parts = line.split(self.delim)
        if self.header == None:
            header_string = self.__generate_header_from_line__(line)
            print(f"Warning: Adding line ({line}) with unset header. Creating a default header ({header_string}).")
            self.header = header_string
        if len(parts) != len(self.header.split(',')):
            return False
        return True
    
    def __verify_header__(self, header):
        """
        Verify that the header is in the correct format.
        :param header: The header to verify.
        :return: True if the header is in the correct format, False otherwise.
        """
        parts = header.split(self.delim)
        if header[0] != "#":
            return False
        if len(parts) < 3:
            raise ValueError(f"Line {header} does not contain enough fields.")
        return True
    
    def __get_delim__(self, line):
        """
        Get the delimiter used in the line.
        :param line: The line to check.
        :return: The delimiter used in the line.
        """
        if len(line.split()) >= 3:
            if len(line.split(',')) >= 3:
                raise ValueError("Line contains both comma and tab delimiters.")
            return None # This acts as whitespace
        elif len(line.split(',')) >= 3:
            return ','
        else:
            raise ValueError("Line does not contain a valid delimiter (comma or tab).")

    def __generate_header_from_line__(self, line):
        """
        Generate a generic header from the line.
        :param line: The line to generate the header from.
        :return: The generated header.
        """
        header_string = "#All_Accessions,Match_Label,Match_E_Value"
        if line is not None:
            parts = line.split(self.delim)
            if len(parts) > 3:
                # File has extra columns
                for i in range(len(parts[3:])):
                    header_string += f",Tier_{i}"
        return header_string
    
    def __len__(self):
        return len(self.sorted_list)

    def output_string(self):
        """
        Output the string representation of the label file.
        :return: The string representation of the label file.
        """
        mystr = self.header + "\n"
        for line in self.sorted_list:
            mystr += line.output_string() + "\n"

        return mystr.rstrip()

    def load_from_file(self, filename):
        """
        Load the label file into the sorted list.
        :param filename: The name of the file to load.
        return: None, but populates the sorted list.
        """
        if not path.exists(filename):
            raise ValueError(f"Label file {filename} does not exist.")
        
        elif path.getsize(filename) > self.max_file_size:
            raise ValueError(f"Label file {filename} is too large ({self.__human_readable_size__(path.getsize(filename))}). Please check the file size. It should be less than {self.__human_readable_size__(self.max_file_size)}.")
        
        elif path.getsize(filename) == 0:
            # File is empty, just set a default header
            self.header = self.__generate_header_from_line__(None)
        
        else:
            with open(filename, "r") as f:
                counter = -1

                for line in f:
                    counter += 1
                    line = line.strip()

                    # Check the first line for the header
                    if counter == 0:
                        self.delim = self.__get_delim__(line)

                        is_header = self.__verify_header__(line)

                        if is_header:
                            # If the header is valid, set it
                            self.header = ','.join(line.split(self.delim))
                        else:
                            # If the header is not valid, build it from the first line, then add the line
                            self.header = self.__generate_header_from_line__(line)

                            # Add the line to the sorted list
                            self.add_line(line)

                    else:
                        new_delim = self.__get_delim__(line)
                        if new_delim != self.delim:
                            raise ValueError(f"Label file {filename} has inconsistent delimiters. Expected {self.delim} but found {new_delim}.")
                    
                        parts = line.split(self.delim)

                        if len(parts) < 3:
                            raise ValueError(f"Label file {filename} is not in the correct format. Line index {counter} has too few items when split on {self.delim if self.delim == ',' else 'tab'}: {line}")
                        
                        # Add the line to the sorted list
                        self.add_line(line)

    def save_to_file(self, filename):
        with open(filename, "w") as f:
            f.write(self.header + "\n")
            for record in self:
                key = record.key
                rest = record.rest
                f.write(f"{key},{rest}\n")

    # Utility methods to get, add, and update records    
    def get_record(self, key):
        """
        Get a record by its key.
        :param key: The key to search for.
        :return: The LabelRecord object if found, None otherwise.
        """
        tmp_record = LabelRecord(key, "")

        # Use bisect_left to find the index of the key
        index = self.sorted_list.bisect_left(tmp_record)
        if index < len(self.sorted_list) and self.sorted_list[index].key == key:
            return self.sorted_list[index]
        return None
    
    def add_line(self, line):
        """
        Add a label file line to the sorted list.
        :param line: The line to add.
        """
        line_delim = self.__get_delim__(line)
        if self.delim != line_delim:
            raise ValueError(f"Label file {self.filename} has inconsistent delimiters. Expected {self.delim} but found {line_delim}.")
            
        # Check if the line is in the correct format
        if not self.__verify_line_vs_header__(line):
            raise ValueError(f"Line '{line}' does not match the header format: '{self.header}'")

        parts = line.split(self.delim)

        key = parts[0]
        rest = ','.join(parts[1:])

        new_record = LabelRecord(key, rest)

        # Check if key is already in the LabelFile
        existing = self.get_record(key)

        if existing is not None:
            # If the key is already in the SortedList, we need to compare the e-values
            if existing.e_value_le(new_record):
                # Old record has a smaller (or equal) e-value, so skip
                return
            else:
                # New record has a smaller e-value, so we need to replace the old one
                # Only replace the start of the rest (leave the tiers if they exist)
                existing_rest = existing.rest.split(',')
                new_rest = new_record.rest.split(',')
                existing.rest = ','.join(new_rest[:2] + existing_rest[2:])

                return

        else:
            self.sorted_list.add(new_record)
            return

    def fetch_records_by_label(self, label, n=0):
        """
        O(N)
        Fetch all records with a specific label.
        :param label: The label to search for.
        :param n: The class column to search for (default=0, for first class).
        :return: A SortedList of LabelRecord objects with the specified label.
        """
        records = SortedList()
        # Record.rest has the format "label,e_value,class_0_label,class_0_colour,class_1_label,class_1_colour,..."
        idx = n # index of the class label we're interested in

        for record in self:
            split = record.rest.split(',')
            if split[idx] == label:
                records.add(record)

        return records

    def add_class_column(self, *names):
        """
        O(N)
        Add new columns to the label file.
        :param names: The names of the new columns to add. (added to the header)
        :return: None
        """
        if self.header is None:
            raise ValueError("Header is not set. Cannot add class columns.")

        # Add the new column names to the header
        for name in names:
            self.header += "," + name

        # Add the new column names to each record in the sorted list
        for record in self:
            record.rest += "," + ",".join([""] * len(names))

    def remove_class_column(self, name):
        """
        O(N)
        Remove a column from the label file.
        :param name: The name of the column (in the header) to remove.
        :return: True if successful, False otherwise.
        """

        if self.header is None:
            raise ValueError("Header is not set. Cannot remove class columns.")

        # Remove the column name from the header
        parts = self.header.split(',')
        idx = 3
        while idx < len(parts):
            if parts[idx] == name:
                break
            idx += 1

        if idx >= len(parts):
            return False
        
        # Remove the column name from the header
        del parts[idx]

        # Update header
        self.header = ','.join(parts)

        # For each record, remove the column
        for record in self:
            split = record.rest.split(',')
            del split[idx]
            record.rest = ','.join(split)

        return True
    

if __name__ == "__main__":
    print("Beginning LabelFile tests...")

    passed = 0
    failed = 0

    # 1 Check the max file size
    filename = "lf_testing/l1_new.label"
    try:
        label_file = LabelFile(filename, max_file_size=1)
    except ValueError as e:
        if str(e).startswith(f"Label file {filename} is too large"):
            passed += 1
        else:
            failed += 1
            print(f"Failed to check max file size: {e}")

    # 2 Check invalid path
    filename = "lf_testing/NON_EXSISTENT.label.blah.foo"
    try:
        label_file = LabelFile(filename)
    except ValueError as e:
        if str(e).startswith(f"Label file {filename} does not exist"):
            passed += 1
        else:
            failed += 1
            print(f"Failed to check invalid path: {e}")

    # 3 Check add_line and __str__
    line_3 = "accession1,label1,0.001"
    label_file = LabelFile()
    try:
        label_file.add_line(line_3)
        if str(label_file) == f"LabelFile(None)\n({line_3})":
            passed += 1
        else:
            failed += 1
            print(f"Failed to add line: {line_3}.\nExpected:\nLabelFile(None)\n({line_3})\nGot:\n{str(label_file)}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #3: {e}")

    # 4 Add another line
    line_4 = "accession3,label2,0.002"
    try:
        label_file.add_line(line_4)
        if str(label_file) == f"LabelFile(None)\n({line_3})\n({line_4})":
            passed += 1
        else:
            failed += 1
            print(f"Failed to add line: {line_4}.\nExpected:\nLabelFile(None)\n({line_3})\n({line_4})\nGot:\n{str(label_file)}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #4: {e}")

    # 5 Add another line between the two
    line_5 = "accession2,label3,0.001"
    try:
        label_file.add_line(line_5)
        if str(label_file) == f"LabelFile(None)\n({line_3})\n({line_5})\n({line_4})":
            passed += 1
        else:
            failed += 1
            print(f"Failed to add line: {line_5}.\nExpected:\nLabelFile(None)\n({line_3})\n({line_5})\n({line_4})\nGot:\n{str(label_file)}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #5: {e}")
    
    # 6 Try add a line with a larger e-value
    line_6 = "accession1,label4,0.002"
    try:
        label_file.add_line(line_6)
        if str(label_file) == f"LabelFile(None)\n({line_3})\n({line_5})\n({line_4})":
            passed += 1
        else:
            failed += 1
            print(f"Failed to add line: {line_6}.\nExpected:\nLabelFile(None)\n({line_3})\n({line_5})\n({line_4})\nGot:\n{str(label_file)}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #6: {e}")

    # 7 Add a line with a smaller e-value
    line_7 = "accession1,label5,0.00001"
    try:
        label_file.add_line(line_7)
        if str(label_file) == f"LabelFile(None)\n({line_7})\n({line_5})\n({line_4})":
            passed += 1
        else:
            failed += 1
            print(f"Failed to add line: {line_7}.\nExpected:\nLabelFile(None)\n({line_7})\n({line_5})\n({line_4})\nGot:\n{str(label_file)}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #7: {e}")

    # 8 Check get_record with a valid key (middle)
    key_record = LabelRecord("accession2", "")
    try:
        record = label_file.get_record(key_record.key)
        if record is not None and record.key == key_record.key:
            passed += 1
        else:
            failed += 1
            print(f"Failed to get record with key: {key_record.key}.\nExpected:\n{key_record}\nGot:\n{record}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #8: {e}")

    # 9 Check get_record with a valid key (first)
    key_record = LabelRecord("accession1", "")
    try:
        record = label_file.get_record(key_record.key)
        if record is not None and record.key == key_record.key:
            passed += 1
        else:
            failed += 1
            print(f"Failed to get record with key: {key_record.key}.\nExpected:\n{key_record}\nGot:\n{record}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #9: {e}")

    # 10 Check get_record with a valid key (last)
    key_record = LabelRecord("accession3", "")
    try:
        record = label_file.get_record(key_record.key)
        if record is not None and record.key == key_record.key:
            passed += 1
        else:
            failed += 1
            print(f"Failed to get record with key: {key_record.key}.\nExpected:\n{key_record}\nGot:\n{record}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #10: {e}")

    # 11 Check get_record with an invalid key
    key_record = LabelRecord("accession4", "")
    try:
        record = label_file.get_record(key_record.key)
        if record is None:
            passed += 1
        else:
            failed += 1
            print(f"Failed to get record with key: {key_record.key}.\nExpected:\nNone\nGot:\n{record}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #11: {e}")

    # 12 Read a label file without a header in the new format (1 additional columns)
    filename = "lf_testing/lf_12.csv"
    expected = ''
    filecontents = ''
    with open (filename, "r") as f:
        for line in f:
            expected += "(" + line.strip() + ")\n"
            filecontents += line.strip() + "\n"
    expected = expected.strip()
    filecontents = filecontents.strip()
    expected = f"LabelFile({filename})\n{expected}"

    try:
        label_file = LabelFile(filename)
        if str(label_file) == expected:
            passed += 1
        else:
            failed += 1
            print(f"Failed to read label file without a header in the new format.\nExpected:\n{expected}\nGot:\n{str(label_file)}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #12: {e}")
        
    # 13 Check (12) against the outstring
    try:
        label_file = LabelFile(filename)
        expected = f"#Accession,Label,E_Value,Class_0,Class_0_Colour\n{filecontents}"
        if label_file.output_string() == expected:
            passed += 1
        else:
            failed += 1
            print(f"Failed to check output string against expected.\nExpected:\n{expected}\nGot:\n{label_file.output_string()}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #13: {e}")

    # 14 Read a label file without a header in the new format (2 additional columns)
    filename = "lf_testing/lf_14.csv"
    filecontents = ''
    with open (filename, "r") as f:
        for line in f:
            filecontents += line.strip() + "\n"
    filecontents = filecontents.strip()

    try:
        label_file = LabelFile(filename)
        expected = f"#Accession,Label,E_Value,Class_0,Class_0_Colour,Class_1,Class_1_Colour\n{filecontents}"
        if label_file.output_string() == expected:
            passed += 1
        else:
            failed += 1
            print(f"Failed to read label file without a header (2 addtl columns) in the new format.\nExpected:\n{expected}\nGot:\n{label_file.output_string()}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #14: {e}")

    # 15 Write (14) to a file
    outfilename = "lf_testing/lf_15.csv"
    infilename = "lf_testing/lf_14.csv"
    try:
        label_file = LabelFile(infilename)
        label_file.save_to_file(outfilename)
        with open(outfilename, "r") as f1:
            with open(infilename, "r") as f2:
                lines1 = f1.readlines()
                lines2 = f2.readlines()
                lines1 = [line.strip() for line in lines1]
                lines2 = [line.strip() for line in lines2]

                lines2.insert(0, "#Accession,Label,E_Value,Class_0,Class_0_Colour,Class_1,Class_1_Colour")
                
                if lines1 == lines2:
                    passed += 1
                else:
                    failed += 1
                    print(f"Failed to save label file to {outfilename}.\nExpected:\n{lines2}\nGot:\n{lines1}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #15: {e}")

    # 16 Read a label file without a header in the new format (0 additional columns)
    filename = "lf_testing/lf_16.csv"
    expected = ''
    with open(filename, "r") as f:
        for line in f:
            expected += line.strip() + "\n"
    expected = expected.strip()
    expected = f"#Accession,Label,E_Value\n{expected}"

    try:
        label_file = LabelFile(filename)
        if label_file.output_string() == expected:
            passed += 1
        else:
            failed += 1
            print(f"Failed to read label file without a header in the new format.\nExpected:\n{expected}\nGot:\n{label_file.output_string()}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #16: {e}")

    # 17 Read a label file with a header in the new format (0 additional columns)
    filename = "lf_testing/lf_17.csv"
    expected = ''
    with open(filename, "r") as f:
        expected = f.read()
    expected = expected.strip()

    try:
        label_file = LabelFile(filename)
        if label_file.output_string() == expected:
            passed += 1
        else:
            failed += 1
            print(f"Failed to read label file with a header in the new format.\nExpected:\n{expected}\nGot:\n{label_file.output_string()}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #17: {e}")

    # 18 Read a label file with a header in the old format (csv header?) - Inconsistent delimiters!
    filename = "lf_testing/lf_18.csv"

    try:
        label_file = LabelFile(filename)
        failed += 1
        print(f"No exception raised on test #18. Expected a ValueError due to inconsistent delimiters.")
    except ValueError as e:
        if str(e).startswith(f"Label file {filename} has inconsistent delimiters. Expected"):
            passed += 1
        else:
            failed += 1
            print(f"Failed to check inconsistent delimiters on test #18.\nExpected:\nValueError\nGot:\n{e}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #18: {e}")

    # 19 Read a label file with a header in the old format (tab header?)
    filename = "lf_testing/lf_19.tsv"
    expected = ''
    with open(filename, "r") as f:
        for line in f:
            expected += ','.join(line.strip().split()) + "\n"
    expected = expected.strip()

    try:
        label_file = LabelFile(filename)
        if label_file.output_string() == expected:
            passed += 1
        else:
            failed += 1
            print(f"Failed to read label file with a header in the old format.\nExpected:\n{expected}\nGot:\n{label_file.output_string()}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #19: {e}")

    # 20 Read a label file without a header in the old format
    filename = "lf_testing/lf_20.tsv"
    expected = ''
    with open(filename, "r") as f:
        for line in f:
            expected += ','.join(line.strip().split()) + "\n"
    expected = expected.strip()
    expected = f"#Accession,Label,E_Value\n{expected}"

    try:
        label_file = LabelFile(filename)
        if label_file.output_string() == expected:
            passed += 1
        else:
            failed += 1
            print(f"Failed to read label file without a header in the old format.\nExpected:\n{expected}\nGot:\n{label_file.output_string()}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #20: {e}")

    # 21 Fetch all records with class label 'label1' in column 0
    filename = "lf_testing/lf_21.csv"

    try:
        label_file = LabelFile(filename)
        records = label_file.fetch_records_by_label("class1_label_1", 0)
        if len(records) == 2:
            passed += 1
        else:
            failed += 1
            print(f"Failed to fetch records with class label 'label1' in column 0.\nExpected:\n2 records\nGot:\n{len(records)}")
    except Exception as e:
        failed += 1
        print(f"Unexpected exception on test #21: {e}")

    print(f"Passed {passed} out of {passed + failed} tests ({round(passed/(passed+failed)*100, 2)}%).")
