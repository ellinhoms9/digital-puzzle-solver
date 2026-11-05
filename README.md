# Digital Puzzle Solver

## Description

This Python program finds the longest possible sequence of numeric fragments.

Two fragments can be connected when the last two digits of the first fragment match the first two digits of the second fragment.

Each fragment can be used only once.

The program uses mathematical optimization to maximize the number of connected fragments.

## Requirements

- Python 3
- NumPy
- SciPy

## Installation

1. Install Python 3.
2. Download or clone this project.
3. Open a terminal in the project directory.
4. Install the required libraries:

```bash
python3 -m pip install -r requirements.txt
```

## Input file

Place the numeric fragments in `source.txt`.

Each fragment must be on a separate line and contain only digits.

Example:

```text
608017
248460
962282
994725
177092
```

You can replace the contents of `source.txt` with another dataset.

## Running the program

Run the following command:

```bash
python3 main.py
```

The program displays:
- The total number of input fragments.
- The number of fragments used in the solution.
- The length of the resulting sequence.
- The resulting numeric sequence.
- The order of the connected fragments.

## Dependencies and licensing

The program uses Python, NumPy and SciPy. These are legally available software tools. NumPy and SciPy are open-source projects.

No proprietary or pirated software is required.