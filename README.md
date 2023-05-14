# DLX Sudoku Solver

## Table of Contents

- [Description](#description)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## Description

This is a C++ implementation of Donald Knuth's Algorithm DLX which enumerates all solutions to an exact cover problem. It is used to solve Sudoku puzzles.

## Installation

Change the file locations in the global variable area (around line 215) to something that works for you and compile the source then put the input file (puzzles.txt) and the executable in your storage device. 

## Usage

This is a CLI program with no input parameters. The original and any solutions go to STDOUT, Error messages go to STDERR. Be advised that lines 13 and 14 of puzzles will cause the program to search for solutions until the cows come home. I included them in case somebody is bored. ;<) 

## Contributing

Send me a pull request. I ask that you please reach out first. Thank you.

## License

N/A

## Acknowledgements
* Dr Knuth's *Dancing Links* Paper: http://www.ocf.berkeley.edu/~jchu/publicportal/sudoku/0011047.pdf
* Jonathan Chu's detailed explanation on turning a Sudoku grid into an exact cover problem: https://www.ocf.berkeley.edu/~jchu/publicportal/sudoku/sudoku.paper.html
