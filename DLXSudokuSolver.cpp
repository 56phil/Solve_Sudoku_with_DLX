/*
 A slightly safer sudoku solver

 This program is a modification of Karl Hajal's work done in 2017. The structure
 named Node was converted to a class with accessors and mutators and the name
 was changed to DLX_Node] as suggested by ChatGPT. The actual solver that Karl
 did (logic and code organization) is mostly untouched.

 ===============================================================================
 ChatGPT summary follows...
 ===============================================================================

The code consists of several parts:

1. Definitions and global variables: The necessary C++ libraries are included,
and constants and global variables are defined.

2. The `DLX_Node` class: This class represents a node in the Dancing Links data
structure. It has member variables and accessor/mutator methods to manipulate
the node's properties.

3. `coverColumn` and `uncoverColumn` functions: These functions are used to
cover and uncover columns in the Dancing Links matrix.

4. `findSolution` function: This is the main backtracking algorithm that finds
all solutions to the Sudoku puzzle. It chooses a column to cover
deterministically based on the size of the column. It recursively explores all
possible combinations of rows until a valid solution is found.

5. `BuildSparseMatrix` function: This function builds the initial sparse matrix
representation of the Sudoku puzzle. It sets up the constraints for the exact
cover problem.

6. `BuildLinkedList` function: This function converts the sparse matrix
representation into a toroidal doubly linked list, which is the data structure
used in the Dancing Links algorithm. Each node in the list represents a 1 in the
sparse matrix.

7. `TransformListToCurrentGrid` function: This function covers the nodes in the
list that correspond to the values already present in the Sudoku grid. This step
eliminates the possibilities for those cells and reduces the search space.

8. Helper functions: There are several helper functions for printing the Sudoku
grid and converting between matrix and string representations of the puzzle.

9. The `main` function: The program reads Sudoku puzzles from a file or uses a
default puzzle if the file cannot be opened. It then solves each puzzle using
the `solvePuzzle` function. The execution time for each puzzle is measured using
the `std::chrono` library. Finally, the solutions are written to files.

The Dancing Links algorithm is an efficient algorithm for solving exact cover
problems, and it has been adapted here to solve Sudoku puzzles. The algorithm
recursively explores the solution space and finds all possible exact covers,
which correspond to valid solutions to the Sudoku puzzle./ This program is a
slightly modified version of a Sudoku solver. The solver is based on the Dancing
Links algorithm, also known as Algorithm X, which is an efficient algorithm for
solving exact cover problems. The main idea behind the algorithm is to represent
the Sudoku puzzle as an exact cover problem and then use a backtracking
algorithm to find the solutions.

 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ChatGPT's exact cover problem defined:
 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 An exact cover problem is a type of combinatorial optimization problem that
 involves finding a subset of a given collection of sets, called the "covering"
 sets, such that each element appears exactly once in the selected sets. In
 other words, the goal is to find a collection of sets that covers every element
 exactly once.

 Formally, an exact cover problem can be defined as follows: Given a universe U
 and a collection S of subsets of U, the problem is to find a subcollection S'
 of S such that every element in U appears in exactly one set in S'. The
 subcollection S' is an exact cover if and only if each element in U appears in
 exactly one set in S'.

 Exact cover problems can be represented as an exact cover matrix, where each
 row represents an element from the universe U and each column represents a
 subset from the collection S. The matrix contains binary values indicating
 whether an element is present in a particular subset.

 Solving an exact cover problem involves finding a solution that satisfies the
 exact cover condition. This can be achieved through various algorithms, such as
 the Dancing Links algorithm, which efficiently explores the solution space and
 finds all possible exact covers.

 Exact cover problems have applications in various areas, including computer
 science, mathematics, and operations research. They are used to solve puzzles,
 scheduling problems, packing problems, and other optimization and constraint
 satisfaction problems. One famous example of an exact cover problem is the
 Sudoku puzzle, where the goal is to find an exact cover of a 9x9 grid with
 certain constraints on row, column, and box completeness.
 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <string>
#include <vector>

struct separated : std::numpunct<char> {
  std::string do_grouping() const { return "\03"; }
};

struct puzzleStruct {
  std::string rawPuzzle;
  std::vector<std::string> solutionStrings;
  double timeToSolve;
  int state;
};

enum {
  OK = 0,
  fileOpenError,
  rawPuzzleSizeError,
  mat2strErr,
  str2matErr,
  noSolution,
};

class DLX_Node {

public:
  DLX_Node() {
    left = this;
    right = this;
    up = this;
    down = this;
    head = this;
    size = 0;
    candidate = 0;
    row = 1;
    column = 1;
  }

  DLX_Node *getLeft() { return left; }
  void setLeft(DLX_Node *_l) { left = _l; }

  DLX_Node *getRight() { return right; }
  void setRight(DLX_Node *_r) { right = _r; }

  DLX_Node *getUp() { return up; }
  void setUp(DLX_Node *_u) { up = _u; }

  DLX_Node *getDown() { return down; }
  void setDown(DLX_Node *_d) { down = _d; }

  DLX_Node *getHead() { return head; }
  void setHead(DLX_Node *_h) { head = _h; }

  int getSize() { return size; }
  void setSize(int _s) { size = _s; }
  void decSize() { size--; }
  void incSize() { size++; }

  int getCandidate() { return candidate; }
  void setCandidate(int _c) { candidate = _c; }
  void decCandidate() { candidate--; }
  void incCandidate() { candidate++; }

  int getRow() { return row; }
  void setRow(int _r) { row = _r; }
  void decRow() { row--; }
  void incRow() { row++; }

  int getColumn() { return column; }
  void setColumn(int _c) { column = _c; }
  void decColumn() { column--; }
  void incColumn() { column++; }

private:
  DLX_Node *left;
  DLX_Node *right;
  DLX_Node *up;
  DLX_Node *down;
  DLX_Node *head;

  int size; // used for Column header
  int candidate;
  int row;
  int column;
};

//=====================================================================//
// --- Golbal consts & variables ------------------------------------- //
//=====================================================================//
class GlobalData {
public:
  static bool isSolved;
  static double puzzleTimeLimit;
  static int expectedCSVlength;
  static std::chrono::steady_clock::time_point puzzleStartTime;
  static std::string inFileName;
  static std::string outFileName;
  static std::vector<puzzleStruct> puzzleStructs;
  static puzzleStruct currentPuzzleStruct;
};

const int MAX_K(999);
const int SIZE(9);
const int COL_NB(4 * SIZE * SIZE);
const int ROW_NB(SIZE *SIZE *SIZE);
const int SIZE_SQRT((int)sqrt((double)SIZE));
const int SIZE_SQUARED(SIZE *SIZE);

DLX_Node Head;
DLX_Node *HeadNode(&Head);

DLX_Node *solution[MAX_K];
DLX_Node *orig_values[MAX_K];

//=====================================================================//
// --- initialize members of GlobalData ------------------------------ //
//=====================================================================//
std::chrono::steady_clock::time_point GlobalData::puzzleStartTime =
    std::chrono::steady_clock::now();
std::string GlobalData::inFileName = "";
std::string GlobalData::outFileName = "";
double GlobalData::puzzleTimeLimit = 2.5;
bool GlobalData::isSolved = false;
std::vector<puzzleStruct> GlobalData::puzzleStructs;
puzzleStruct GlobalData::currentPuzzleStruct;
int GlobalData::expectedCSVlength = 163;

//=====================================================================//
// --- Prototypes ---------------------------------------------------- //
//=====================================================================//
bool isSolution(const std::string &s);
bool isTooSlow();
int getFileNamesFromCommandLine(int a, char *b[]);
int getPuzzlesFromStorage();
int matrix2string(std::string &result, int mat[][SIZE]);
int string2matrix(std::string const rawSudoku, int mat[][SIZE]);
int writeResults();
std::string number2words(long num);
void mapSolutionToGrid(int sudoku[][SIZE]);
void primaryLoop(int (*activeMatrix)[9], int &totalSolutionsFound);
void printGrid(int sudoku[][SIZE]);
void solvePuzzle(int sudoku[][SIZE]);
void wrapUp(int totalSolutionsFound);

//=====================================================================//
// --- MAIN ---------------------------------------------------------- //
//=====================================================================//
int main(int argc, char *argv[]) {

  std::locale our_local(std::cout.getloc(), new separated);
  std::cout.imbue(our_local);

  getFileNamesFromCommandLine(argc, argv);

  int totalSolutionsFound(0);
  int returnCode(OK);
  std::string defaultPuzzle("00000000000000308500102000000050700000400010009000"
                            "0000500000073002010000000040009");
  int activeMatrix[SIZE][SIZE];

  GlobalData::puzzleStructs.clear();
  if ((returnCode = getPuzzlesFromStorage()) != OK) {
    GlobalData::puzzleStructs.clear();
  }
  primaryLoop(activeMatrix, totalSolutionsFound);
  wrapUp(totalSolutionsFound);
  return returnCode;
}

//=====================================================================//
// --- main loop ----------------------------------------------------- //
//=====================================================================//
void primaryLoop(int (*activeMatrix)[9], int &totalSolutionsFound) {
  int cntr(1);

  for (puzzleStruct &ps : GlobalData::puzzleStructs) {
    GlobalData::puzzleStartTime = std::chrono::steady_clock::now();
    GlobalData::currentPuzzleStruct = ps;
    GlobalData::currentPuzzleStruct.solutionStrings.clear();

    string2matrix(GlobalData::currentPuzzleStruct.rawPuzzle, activeMatrix);
    printGrid(activeMatrix);
    solvePuzzle(activeMatrix);
    GlobalData::currentPuzzleStruct.state =
        GlobalData::isSolved ? OK : noSolution;
    std::chrono::steady_clock::time_point end =
        std::chrono::steady_clock::now();
    ps.timeToSolve = std::chrono::duration_cast<std::chrono::duration<double>>(
                         end - GlobalData::puzzleStartTime)
                         .count();

    long tempLong(GlobalData::currentPuzzleStruct.solutionStrings.size());
    std::string sol1(number2words(tempLong));
    sol1 += " solution";
    sol1 += (tempLong != 1 ? "s" : "");
    sol1 += " identified.\n";

    // if (tempLong == 0) {
    //   sol1 = "No solutions";
    // } else if (tempLong == 1) {
    //   sol1 = "One solution";
    // } else {
    //   sol1 = number2words(tempLong);
    //   sol1 += " solutions";
    // }
    // sol1 += " identified.\n\n\n";

    std::cout << cntr++ << ". Execution time: " << ps.timeToSolve
              << " seconds\n"
              << sol1;

    totalSolutionsFound += tempLong;
    ps = GlobalData::currentPuzzleStruct;
  }
}

//=====================================================================//
// --- wrao up ------------------------------------------------------- //
//=====================================================================//
void wrapUp(int totalSolutionsFound) {
  writeResults();
  std::cout << "\n\n"
            << number2words(GlobalData::puzzleStructs.size())
            << " puzzles processed.\nA total of "
            << number2words(totalSolutionsFound) << " solutions identified.";
}

//=====================================================================//
// --- DLX Functions ------------------------------------------------- //
//=====================================================================//
void coverColumn(DLX_Node *col) {
  col->getLeft()->setRight(col->getRight());
  col->getRight()->setLeft(col->getLeft());
  for (DLX_Node *node(col->getDown()); node != col; node = node->getDown()) {
    for (DLX_Node *temp(node->getRight()); temp != node;
         temp = temp->getRight()) {
      temp->getDown()->setUp(temp->getUp());
      temp->getUp()->setDown(temp->getDown());
      temp->getHead()->decSize();
    }
  }
}

void uncoverColumn(DLX_Node *col) {
  for (DLX_Node *node(col->getUp()); node != col; node = node->getUp()) {
    for (DLX_Node *temp(node->getLeft()); temp != node;
         temp = temp->getLeft()) {
      temp->getHead()->incSize();
      temp->getDown()->setUp(temp);
      temp->getUp()->setDown(temp);
    }
  }
  col->getLeft()->setRight(col);
  col->getRight()->setLeft(col);
}

void findSolution(int kount) {
  if (isTooSlow()) {
    return;
  }

  if (HeadNode->getRight() == HeadNode) {
    int Grid[SIZE][SIZE] = {{0}};
    mapSolutionToGrid(Grid);
    printGrid(Grid);
    GlobalData::isSolved = true;
    return;
  }

  // Choose Column Object Deterministically: Choose the column with the smallest
  // Size
  DLX_Node *Col(HeadNode->getRight());
  for (DLX_Node *temp(Col->getRight()); temp != HeadNode;
       temp = temp->getRight())
    if (temp->getSize() < Col->getSize())
      Col = temp;

  coverColumn(Col);

  for (DLX_Node *temp(Col->getDown()); temp != Col; temp = temp->getDown()) {
    solution[kount] = temp;
    for (DLX_Node *node(temp->getRight()); node != temp;
         node = node->getRight()) {
      coverColumn(node->getHead());
    }

    findSolution(kount + 1);

    temp = solution[kount];
    solution[kount] = nullptr;
    Col = temp->getHead();
    for (DLX_Node *node(temp->getLeft()); node != temp;
         node = node->getLeft()) {
      uncoverColumn(node->getHead());
    }
  }

  uncoverColumn(Col);
}

//=====================================================================//
// --- Functions to turn a Sudoku grid into an Exact Cover problem --- //
//=====================================================================//
// --- BUILD THE INITIAL MATRIX CONTAINING ALL POSSIBILITIES --------- //
void BuildSparseMatrix(bool matrix[ROW_NB][COL_NB]) {

  // Constraint 1: There can only be one value in any given cell
  int j(0), counter(0);
  for (int i(0); i < ROW_NB; i++) { // iterate over all rows
    matrix[i][j] = 1;
    counter++;
    if (counter >= SIZE) {
      j++;
      counter = 0;
    }
  }

  // Constraint 2: There can only be one instance of a number in any given row
  int x(0);
  counter = 1;
  for (j = SIZE_SQUARED; j < 2 * SIZE_SQUARED; j++) {
    for (int i(x); i < counter * SIZE_SQUARED; i += SIZE)
      matrix[i][j] = 1;

    if ((j + 1) % SIZE == 0) {
      x = counter * SIZE_SQUARED;
      counter++;
    } else
      x++;
  }

  // Constraint 3: There can only be one instance of a number in any given
  // column
  j = 2 * SIZE_SQUARED;
  for (int i(0); i < ROW_NB; i++) {
    matrix[i][j] = 1;
    j++;
    if (j >= 3 * SIZE_SQUARED)
      j = 2 * SIZE_SQUARED;
  }

  // Constraint 4: There can only be one instance of a number in any given
  // region
  x = 0;
  for (j = 3 * SIZE_SQUARED; j < COL_NB; j++) {

    for (int l(0); l < SIZE_SQRT; l++) {
      for (int k(0); k < SIZE_SQRT; k++)
        matrix[x + l * SIZE + k * SIZE_SQUARED][j] = 1;
    }

    int temp = j + 1 - 3 * SIZE_SQUARED;

    if (temp % (int)(SIZE_SQRT * SIZE) == 0)
      x += (SIZE_SQRT - 1) * SIZE_SQUARED + (SIZE_SQRT - 1) * SIZE + 1;
    else if (temp % SIZE == 0)
      x += SIZE * (SIZE_SQRT - 1) + 1;
    else
      x++;
  }
}

//=====================================================================//
// --- BUILD A TOROIDAL DOUBLY LINKED LIST OUT OF THE SPARSE MATRIX ---//
//=====================================================================//
void BuildLinkedList(bool matrix[ROW_NB][COL_NB]) {

  DLX_Node *header(new DLX_Node);
  DLX_Node *temp(header);

  // Create all Column Nodes
  for (int i(0); i < COL_NB; i++) {
    DLX_Node *newNode(new DLX_Node);
    newNode->setRight(header);
    newNode->setLeft(temp);
    temp->setRight(newNode);
    temp = newNode;
  }

  int candidate(0);
  int row(1);
  int column(1);
  // Add a Node for each 1 present in the sparse matrix and update Column Nodes
  // accordingly
  for (int i(0); i < ROW_NB; i++) {
    DLX_Node *top(header->getRight());
    DLX_Node *prev(nullptr);

    if (i != 0 && i % SIZE_SQUARED == 0) {
      candidate -= SIZE - 1;
      row++;
      column -= SIZE - 1;
    } else if (i != 0 && i % SIZE == 0) {
      candidate -= SIZE - 1;
      column++;
    } else {
      candidate++;
    }

    for (int j(0); j < COL_NB; j++, top = top->getRight()) {
      if (matrix[i][j]) {
        DLX_Node *newNode(new DLX_Node);
        newNode->setCandidate(candidate);
        newNode->setRow(row);
        newNode->setColumn(column);
        if (prev == nullptr) {
          prev = newNode;
          prev->setRight(newNode);
        }
        newNode->setLeft(prev);
        newNode->setRight(prev->getRight());
        newNode->getRight()->setLeft(newNode);
        prev->setRight(newNode);
        newNode->setHead(top);
        newNode->setDown(top);
        newNode->setUp(top->getUp());
        top->getUp()->setDown(newNode);
        top->incSize();
        top->setUp(newNode);
        if (top->getDown() == top)
          top->setDown(newNode);
        prev = newNode;
      }
    }
  }
  HeadNode = header;
}

//=====================================================================//
// --- COVERS VALUES THAT ARE ALREADY PRESENT IN THE GRID ------------ //
//=====================================================================//
void TransformListToCurrentGrid(int Puzzle[][SIZE]) {
  int index(0);
  for (int i(0); i < SIZE; i++)
    for (int j(0); j < SIZE; j++)
      if (Puzzle[i][j] > 0) {
        DLX_Node *Col(nullptr);
        DLX_Node *temp(nullptr);
        for (Col = HeadNode->getRight(); Col != HeadNode;
             Col = Col->getRight()) {
          for (temp = Col->getDown(); temp != Col; temp = temp->getDown())
            if (temp->getCandidate() == Puzzle[i][j] &&
                (temp->getRow() - 1) == i && (temp->getColumn() - 1) == j)
              goto ExitLoops;
        }
      ExitLoops:
        coverColumn(Col);
        orig_values[index] = temp;
        index++;
        for (DLX_Node *node(temp->getRight()); node != temp;
             node = node->getRight()) {
          coverColumn(node->getHead());
        }
      }
}

//=====================================================================//
// --- Print Functions ----------------------------------------------- //
//=====================================================================//
void mapSolutionToGrid(int Sudoku[][SIZE]) {

  for (int i(0); solution[i] != nullptr; i++) {
    Sudoku[solution[i]->getRow() - 1][solution[i]->getColumn() - 1] =
        solution[i]->getCandidate();
  }

  for (int i(0); orig_values[i] != nullptr; i++) {
    Sudoku[orig_values[i]->getRow() - 1][orig_values[i]->getColumn() - 1] =
        orig_values[i]->getCandidate();
  }
}

//=====================================================================//
// --- PRINTS A SUDOKU GRID ------------------------------------------ //
//=====================================================================//
void printGrid(int Sudoku[][SIZE]) {
  std::string tempString("");

  if ((matrix2string(tempString, Sudoku)) != OK) {
    std::cerr << "matrix2string error: " << tempString << '\n';
    return;
  }

  if (isSolution(tempString) == true) {
    GlobalData::currentPuzzleStruct.solutionStrings.push_back(tempString);
  }

  std::string ext_border = "+", int_border = "|";
  int counter(1);
  int additional(0);
  if (SIZE > 9)
    additional = SIZE;
  for (int i(0); i < ((SIZE + SIZE_SQRT - 1) * 2 + additional + 1); i++) {
    ext_border += '-';

    if (i > 0 && i % ((SIZE_SQRT * 2 + SIZE_SQRT * (SIZE > 9) + 1) * counter +
                      counter - 1) ==
                     0) {
      int_border += '+';
      counter++;
    } else
      int_border += '-';
  }
  ext_border += '+';
  int_border += "|";

  std::cout << ext_border << std::endl;
  for (int i(0); i < SIZE; i++) {
    std::cout << "| ";
    for (int j(0); j < SIZE; j++) {
      if (Sudoku[i][j] == 0)
        std::cout << ". ";
      else
        std::cout << Sudoku[i][j] << " ";
      if (additional > 0 && Sudoku[i][j] < 10)
        std::cout << " ";
      if ((j + 1) % SIZE_SQRT == 0)
        std::cout << "| ";
    }
    std::cout << std::endl;
    if ((i + 1) % SIZE_SQRT == 0 && (i + 1) < SIZE)
      std::cout << int_border << std::endl;
  }
  std::cout << ext_border << std::endl << std::endl;
}

void solvePuzzle(int Sudoku[][SIZE]) {
  bool matrix[ROW_NB][COL_NB] = {{0}};
  BuildSparseMatrix(matrix);
  BuildLinkedList(matrix);
  TransformListToCurrentGrid(Sudoku);
  findSolution(0);
  // if (!GlobalData::isSolved) {
  //   std::cout << "No Solution!" << std::endl;
  // }
  GlobalData::isSolved = false;
}

int writeResults() {
  std::ofstream outFile;
  outFile.open(GlobalData::outFileName.c_str());
  for (puzzleStruct ps : GlobalData::puzzleStructs) {
    outFile << ps.rawPuzzle << '\n';
    for (std::string solution : ps.solutionStrings) {
      outFile << solution << '\n';
    }
  }
  return OK;
}

int getPuzzlesFromStorage() {
  // one line one puzzle. line size must equal SIZE_SQUARED or 2SIZE_SQUARED
  // + 1.
  std::ifstream inFile;
  inFile.open(GlobalData::inFileName.c_str());
  if (inFile.is_open()) {
    puzzleStruct temp;
    std::string line("");
    while (getline(inFile, line)) {
      GlobalData::currentPuzzleStruct.solutionStrings.clear();
      if (line.size() == GlobalData::expectedCSVlength) { // sudoku,solution
        temp.solutionStrings.push_back(
            line.substr(SIZE_SQUARED + 1, SIZE_SQUARED));
        line.resize(81);
      }
      if (line.size() == SIZE_SQUARED) {
        // rawPuzzles.push_back(line);
        temp.rawPuzzle = line;
        GlobalData::puzzleStructs.emplace_back(temp);
      }
    }
    inFile.close();
  } else {
    return fileOpenError;
  }
  return OK;
}

int string2matrix(std::string const rawPuzzle, int mat[SIZE][SIZE]) {
  int i(0), j(0);
  if (rawPuzzle.size() == SIZE_SQUARED) {
    for (char c : rawPuzzle) {
      if (j >= SIZE) {
        i++;
        j = 0;
      }
      mat[i][j++] = c & 0xf;
    }
  }
  return OK;
}

int matrix2string(std::string &result, int mat[][SIZE]) {
  result.clear();
  for (int i(0); i < SIZE; i++) {
    for (int j(0); j < SIZE; j++) {
      result += std::to_string(mat[i][j]);
    }
  }

  if (result.size() != SIZE_SQUARED) {
    return mat2strErr;
  }
  return OK;
}

bool isSolution(const std::string &str) {
  for (int c : str) {
    if ((c & 0xf) == 0) {
      return false;
    }
  }
  return true;
}

bool isTooSlow() {
  auto currentTime = std::chrono::steady_clock::now();

  double elapsed_time =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          currentTime - GlobalData::puzzleStartTime)
          .count();

  // std::cout << "elapsed time: " << elapsed_time << '\n';
  if (elapsed_time > GlobalData::puzzleTimeLimit) {
    GlobalData::isSolved = false;
    std::cerr << "Time limit exceded.\n";
    return true;
  }
  return false;
}

std::string number2words(long num) {
  // Define the word representation for numbers from 0 to 19
  const std::string lessThan20[] = {
      "Zero",    "One",     "Two",       "Three",    "Four",
      "Five",    "Six",     "Seven",     "Eight",    "Nine",
      "Ten",     "Eleven",  "Twelve",    "Thirteen", "Fourteen",
      "Fifteen", "Sixteen", "Seventeen", "Eighteen", "Nineteen"};

  // Define the word representation for multiples of 10 up to 90
  const std::string tens[] = {"",      "",      "Twenty",  "Thirty", "Forty",
                              "Fifty", "Sixty", "Seventy", "Eighty", "Ninety"};

  if (num < 20) {
    return lessThan20[num];
  } else if (num < 100) {
    return tens[num / 10] + (num % 10 != 0 ? " " + lessThan20[num % 10] : "");
  } else if (num < 1000) {
    return lessThan20[num / 100] + " Hundred" +
           (num % 100 != 0 ? " " + number2words(num % 100) : "");
  } else if (num < 1000000) {
    return number2words(num / 1000) + " Thousand" +
           (num % 1000 != 0 ? " " + number2words(num % 1000) : "");
  } else if (num < 1000000000) {
    return number2words(num / 1000000) + " Million" +
           (num % 1000000 != 0 ? " " + number2words(num % 1000000) : "");
  } else {
    return number2words(num / 1000000000) + " Billion" +
           (num % 1000000000 != 0 ? " " + number2words(num % 1000000000) : "");
  }
}

int getFileNamesFromCommandLine(int argc, char *argv[]) {
  // Check for specific parameters
  for (int i(1); i < argc; i++) {
    std::string arg(argv[i]);

    if (arg == "-h" || arg == "--help") {
      // Display help message and exit
      std::cout << "Assume you namrd the program solve...";
      std::cout << "Usage: ./solve <input_file> <output_file>\n";
      return 0;
    } else {
      // Assume it's a positional argument
      std::cout << "Positional argument: " << arg << "\n";
      switch (i) {

      case 1:
        GlobalData::inFileName = arg;
        break;

      case 2:
        GlobalData::outFileName = arg;
        break;
      }
    }
  }

  return OK;
}
// end of source code
