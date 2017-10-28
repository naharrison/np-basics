This is based on an example from the book Clean Code by Robert C. Martin.
The point is that long switch statements (and if/else chains) should never be repeated and should be "used to create polymorphic objects and hidden behind an inheritance relationship so that the rest of the system can't see them."
