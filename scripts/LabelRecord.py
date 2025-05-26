from decimal import Decimal

class LabelRecord:
    """
    Class to store data from each line of a label file efficiently.
    The first element is the key (match accession), and the rest is a comma-separated string of values (rest of label line).
    """
    __slots__ = ('key', 'rest')

    def __init__(self, key, rest):
        self.key = key
        self.rest = rest

    def __str__(self):
        return f"({self.key},{self.rest})"
    
    def output_string(self):
        """
        Returns the raw line that created the LabelRecord, as a String
        """
        return f"{self.key},{self.rest}"
    
    # ----- Methods for SortedList comparison -----
    def __eq__(self, other):
        if isinstance(other, LabelRecord):
            return self.key == other.key and self.rest == other.rest
        return False
    
    def __gt__(self, other):
        if isinstance(other, LabelRecord):
            return self.key > other.key
        return False
    
    def __lt__(self, other):
        if isinstance(other, LabelRecord):
            return self.key < other.key
        return False
    
    def __ge__(self, other):
        if isinstance(other, LabelRecord):
            return not self.key < other.key
        return False
    
    def __le__(self, other):
        if isinstance(other, LabelRecord):
            return not self.key > other.key
        return False
    
    # ----- Checking evalues -----
    def get_e_value(self):
        if self.rest.split(',')[1].strip() == '':
            return Decimal("Infinity")
        return Decimal(self.rest.split(',')[1])

    def e_value_lt(self, other):
        if isinstance(other, LabelRecord):
            return self.get_e_value() < other.get_e_value()
        return False
    
    def e_value_gt(self, other):
        if isinstance(other, LabelRecord):
            return self.get_e_value() > other.get_e_value()
        return False
    
    def e_value_eq(self, other):
        if isinstance(other, LabelRecord):
            return self.get_e_value() == other.get_e_value()
        return False
    
    def e_value_ge(self, other):
        if isinstance(other, LabelRecord):
            return not self.get_e_value() < other.get_e_value()
        return False
    
    def e_value_le(self, other):
        if isinstance(other, LabelRecord):
            return not self.get_e_value() > other.get_e_value()
        return False
    
    # ----- Getting Class Info -----
    def has_classes(self):
        return len(self.rest.split(',')) > 2
    
    def get_num_classes(self):
        return (len(self.rest.split(',')) - 2)
    
    def get_class_n_label(self, n):
        if not self.has_classes():
            return None
        classes = self.rest.split(',')[2:]
        return classes[n]

    def get_class_n_colour(self, n):
        if not self.has_classes():
            return None
        classes = self.rest.split(',')[2:]
        return classes[n * 2 + 1]
    
    # ----- Modifying Class Info -----
    def set_class_n_label(self, n, label):
        """
            n: index of class to set (starts at 0)
        """
        if not self.has_classes():
            return False
        classes = self.rest.split(',')
        classes[n] = label
        self.rest = ','.join(classes)
        return True
    
    # def set_class_n_colour(self, n, colour):
    #     if not self.has_classes() or not self.is_hex_string(colour):
    #         return False
    #     classes = self.rest.split(',')
    #     classes[n * 2 + 1] = colour
    #     self.rest = ','.join(classes)
    #     return True
    
    def pad_new_classes(self, n):
        """
        n: new total number of classes
        """
        classes = self.rest.split(',')
        start, classes = classes[:2], classes[2:]
        
        while len(classes) < n:
            classes += ['']

        if len(classes) > n:
            return False

        self.rest = ','.join(start + classes)
        return True

    # ----- Helper Methods -----
    def is_hex_string(self, s):
        if len(s) != 7:
            return False
        if s[0] != '#':
            return False
        for c in s[1:]:
            if c not in "0123456789abcdefABCDEF":
                return False
        return True