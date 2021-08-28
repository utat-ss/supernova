# A Note on Dynamic 2D Array Formatting

Source: https://stackoverflow.com/questions/40847172/best-way-to-allocate-memory-to-a-two-dimensional-array-in-c

2D arrays are most efficiently stored as a contiguous block, i.e. as a single linear array.

Therefore, the array
```c
[1, 2, 3]
[4, 5, 6]
[7, 8, 9]
```
has memory allocated like the following:
`[1, 2, 3, 4, 5, 6, 7, 8, 9]`

## Array Declaration
Arrays of doubles should be declared like this:
```c
double (*name)[COLUMN_SIZE];
```
Wherein `COLUMN_SIZE` is some constant denoting the column length of the array. For Supernova, we use a column length of 6 (x/y/z/vx/vy/vz)

This way of formatting the array makes it so that the first index of the array yields chunks spaced out by `COLUMN_SIZE`.

## Array Memory Allocation
Memory allocation can be done like the following:
```c
name = malloc(n * sizeof(double[COLUMN_SIZE]));
```
Wherein `n` is the number of rows in the array. Similar procedures are done for reallocating memory, if you need to make the array bigger.

## Accessing Array indices
When the array is initialized like this, you can access it just as you would expect with Python.
```c
printf("%f", name[1][2]);
>>> 6