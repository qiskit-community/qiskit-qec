"""Object for indexing hashable objects."""


class Indexer:
    """Object for indexing hashable objects."""

    def __init__(self):
        """Create an indexer."""
        self.forward_lookup = {}
        self.reverse_lookup = {}
        self.current_index = 0

    def __getitem__(self, key):
        """Forward lookup with automatic insertion of new keys."""
        if key not in self.forward_lookup:
            self.forward_lookup[key] = self.current_index
            self.reverse_lookup[self.current_index] = key
            self.current_index += 1
        return self.forward_lookup[key]

    def rlookup(self, key: int):
        """Reverse lookup."""
        if key not in self.reverse_lookup:
            raise Exception(f"{key} not found")
        return self.reverse_lookup[key]

    def __str__(self):
        """Return string representation."""
        return str(self.forward_lookup)
