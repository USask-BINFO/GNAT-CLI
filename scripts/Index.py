"""
Class to manage the index search for a label file.
"""

from bisect import bisect_left
from os import path

class LabelIndex:
    def __init__(self, label_filename):
        self.accessions = []
        self.offsets = []
        self.label_filename = label_filename
        self.index_filename = label_filename + ".idx"

        if not path.exists(self.index_filename):
            print(f"Index file {self.index_filename} does not exist. Performing regular search.")
            self.index_filename = None
        else:
            self._load_index_data()

    def _load_index_data(self):
        """
        Load the index data from the index file.
        """
        with open(self.index_filename, "r") as index_file:
            # Skip the header lines
            index_file.readline()

            for line in index_file:
                acc, offset = line.strip().split()
                self.accessions.append(acc)
                self.offsets.append(int(offset))

    def _find_offsets(self, target_acc):
        """
        Perform binary search to find the index before the target accession.
        """
        index = bisect_left(self.accessions, target_acc)
        if index == 0:
            # Either accession is the first accession, or not in the file
            start = self.offsets[0]
            end = self.offsets[0]
        elif index == len(self.accessions):
            # Accession is after the last accession in the index file, or not in the file
            start = self.offsets[-1]
            end = None
        else:
            # Accession is between two accessions in the index file
            start = self.offsets[index - 1]
            end = self.offsets[index]

        return start, end
    
    def _find_target_line(self, target, start=None, end=None):
        """
        Find the target line in the label file.
        """
        if self.index_filename is None or start is None:
            return self._bad_search(target) # Fallback to linear search
        
        with open(self.label_filename, 'rb') as label_file:
            label_file.seek(start)

            while True:
                # TODO: Do this in log(n) time
                line = label_file.readline()

                if not line:
                    break # EOF
                
                # Check if we're on the right line
                line_str = line.decode().strip()
                if line_str.split(",")[0] == target:
                    return line_str
                
                pos = label_file.tell()
                if end is not None and pos > end:
                    break # End of the search range
        return None

    def _bad_search(self, target):
        """
        Fallback (linear) search if index file is not available.
        """
        with open(self.label_filename, "r") as label_file:
            label_file.readline()  # Skip header

            for line in label_file:
                line_str = line.strip()
                if line_str.split(",")[0] == target:
                    return line_str
        return None
    
    def search(self, target):
        """
        Search for the target accession in the label file.
        """
        if self.index_filename is None:
            return self._bad_search(target)
        
        start, end = self._find_offsets(target)

        return self._find_target_line(target, start, end)

        