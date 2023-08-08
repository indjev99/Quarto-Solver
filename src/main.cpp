#include <iostream>
#include <cassert>
#include <string>
#include <cctype>
#include <vector>
#include <algorithm>

#define FOR_PROPS(i) for (uint16_t i = 0; i < NUM_PROPS; ++i)
#define FOR_PROPS_VARS(i, j) for (uint16_t i = 0; i < NUM_PROPS; ++i) for (uint16_t j = 0; j < NUM_VARS; ++j) 
#define FOR_PIECES(i) for (uint16_t i = 0; i < NUM_PIECES; ++i)
#define FOR_CELLS(i) for (uint16_t i = 0; i < NUM_CELLS; ++i)
#define FOR_WIN_LEN(i) for (uint16_t i = 0; i < WIN_LEN; ++i)

constexpr uint16_t NUM_VARS = 2u;
constexpr uint16_t NUM_PROPS = 4u;
constexpr uint16_t NUM_PIECES = 1 << NUM_PROPS;
constexpr uint16_t NO_PIECE = NUM_PIECES;

constexpr uint16_t NUM_ROWS = 4;
constexpr uint16_t NUM_COLS = 4;
constexpr uint16_t NUM_CELLS = NUM_ROWS * NUM_COLS;

constexpr uint16_t WIN_LEN = 4;
constexpr uint16_t WIN_SQ_SIDE = 2;

static_assert(WIN_SQ_SIDE * WIN_SQ_SIDE == WIN_LEN);

uint16_t rowColToCell(uint16_t row, uint16_t col)
{
    return row * NUM_COLS + col;
}

uint16_t cellToRow(uint16_t cell)
{
    return cell / NUM_ROWS;
}

uint16_t cellToCol(uint16_t cell)
{
    return cell % NUM_COLS;
}

uint16_t getBit(uint16_t val, uint16_t n)
{
    return (val & (1 << n)) >> n;
}

void setBit(uint16_t& val, uint16_t n)
{
    assert(!getBit(val, n));
    val |= 1 << n;
}

void clearBit(uint16_t& val, uint16_t n)
{
    assert(getBit(val, n));
    val &= ~(1 << n);
}

std::vector<uint16_t> computeWinMasks()
{
    std::vector<uint16_t> winMasks;

    FOR_CELLS(i)
    {
        uint16_t row = cellToRow(i);
        uint16_t col = cellToCol(i);

        if (col + WIN_LEN <= NUM_COLS)
        {
            uint16_t winMask = 0;
            FOR_WIN_LEN(j)
            {
                setBit(winMask, rowColToCell(row, col + j));
            }
            winMasks.push_back(winMask);
        }

        if (row + WIN_LEN <= NUM_ROWS)
        {
            uint16_t winMask = 0;
            FOR_WIN_LEN(j)
            {
                setBit(winMask, rowColToCell(row + j, col));
            }
            winMasks.push_back(winMask);
        }

        if (row + WIN_LEN <= NUM_ROWS && col + WIN_LEN <= NUM_COLS)
        {
            uint16_t winMask = 0;
            FOR_WIN_LEN(j)
            {
                setBit(winMask, rowColToCell(row + j, col + j));
            }
            winMasks.push_back(winMask);
        }

        if (row + WIN_LEN <= NUM_ROWS && col >= WIN_LEN - 1)
        {
            uint16_t winMask = 0;
            FOR_WIN_LEN(j)
            {
                setBit(winMask, rowColToCell(row + j, col - j));
            }
            winMasks.push_back(winMask);
        }

        if (row + WIN_SQ_SIDE <= NUM_ROWS && col + WIN_SQ_SIDE <= NUM_COLS)
        {
            uint16_t winMask = 0;
            FOR_WIN_LEN(j)
            {
                setBit(winMask, rowColToCell(row + j / WIN_SQ_SIDE, col + j % WIN_SQ_SIDE));
            }
            winMasks.push_back(winMask);
        }
    }

    return winMasks;
}

const std::vector<uint16_t> winMasks = computeWinMasks();

struct State
{
    uint16_t movesLeft;
    uint16_t currPiece;
    uint16_t piecesTaken;
    uint16_t cellsTaken;
    uint16_t cellsProps[NUM_PROPS][NUM_VARS];

    State()
    {
        movesLeft = 2 * std::min(NUM_PIECES, NUM_CELLS);
        currPiece = NO_PIECE;
        piecesTaken = 0;
        cellsTaken = 0;

        FOR_PROPS_VARS(i, j)
        {
            cellsProps[i][j] = 0;
        }
    }

    void moveSelect(uint16_t piece)
    {
        assert(!isToPlace());

        setBit(piecesTaken, piece);
        currPiece = piece;

        --movesLeft;
    }

    void undoSelect()
    {
        assert(isToPlace());

        clearBit(piecesTaken, currPiece);
        currPiece = NO_PIECE;

        ++movesLeft;
    }

    void movePlace(uint16_t cell)
    {
        assert(isToPlace());

        setBit(cellsTaken, cell);

        FOR_PROPS(i)
        {
            setBit(cellsProps[i][getBit(currPiece, i)], cell);
        }

        currPiece = NO_PIECE;

        --movesLeft;
    }

    void undoPlace(uint16_t piece, uint16_t cell)
    {
        assert(!isToPlace());

        clearBit(cellsTaken, cell);

        FOR_PROPS(i)
        {
            clearBit(cellsProps[i][getBit(piece, i)], cell);
        }

        currPiece = piece;

        ++movesLeft;
    }

    bool isToPlace()
    {
        return currPiece != NO_PIECE;
    }

    bool isPieceFree(uint16_t piece)
    {
        return !getBit(piecesTaken, piece);
    }

    bool isCellFree(uint16_t cell)
    {
        return !getBit(cellsTaken, cell);
    }

    bool isWon()
    {
        for (uint16_t winMask : winMasks)
        {
            FOR_PROPS_VARS(i, j)
            {
                if ((cellsProps[i][j] & winMask) == winMask) return true;
            }
        }

        return false;
    }

    bool isDone()
    {
        return movesLeft == 0;
    }

    uint16_t getPiece(uint16_t cell)
    {
        if (!getBit(cellsTaken, cell)) return NO_PIECE;
    
        uint16_t piece = 0;
        FOR_PROPS(i)
        {
            if (getBit(cellsProps[i][1], cell)) setBit(piece, i);
        }

        return piece;
    }
};

int16_t INF = 1000;

int16_t sign(int16_t val)
{
    return (0 < val) - (val < 0);
}

int16_t evalSelect(State& state, int16_t alpha, int16_t beta);

int16_t evalPlace(State& state, int16_t alpha, int16_t beta)
{
    int16_t val = -INF;

    uint16_t piece = state.currPiece;

    FOR_CELLS(i)
    {
        if (!state.isCellFree(i)) continue;

        state.movePlace(i);
        int16_t nextVal = evalSelect(state, alpha, beta);
        state.undoPlace(piece, i);

        val = std::max(val, nextVal);
        alpha = std::max(alpha, val);

        if (alpha >= beta) return alpha;
    }

    return val;
}

int16_t evalSelect(State& state, int16_t alpha, int16_t beta)
{
    if (state.isWon()) return state.movesLeft;
    if (state.isDone()) return 0;

    int16_t val = -INF;

    FOR_PIECES(i)
    {
        if (!state.isPieceFree(i)) continue;

        state.moveSelect(i);
        int16_t nextVal = -evalPlace(state, -beta, -alpha);
        state.undoSelect();

        val = std::max(val, nextVal);
        alpha = std::max(alpha, val);

        if (alpha >= beta) return alpha;
    }

    return val;
}

std::string eval(State state)
{
    int16_t val;
    if (state.isToPlace()) val = evalPlace(state, -INF, INF);
    else val = evalSelect(state, -INF, INF);

    if (val == 0) return "Draw";

    std::string str = val > 0 ? "Win in " : "Loss in ";

    str += std::to_string(state.movesLeft - std::abs(val));

    return str;
}

std::string pieceToString(uint16_t piece)
{
    if (piece == NO_PIECE) return "  ";

    char first = getBit(piece, 0) ? 'b' : 'a';
    char second = getBit(piece, 1) ? 'x' : 'o';
    if (getBit(piece, 2)) first = std::toupper(first);
    if (getBit(piece, 3)) second  = std::toupper(second );

    std::string str;
    str += first;
    str += second;

    return str;
}

uint16_t stringToPiece(std::string str)
{
    assert(str.size() == 2);

    if (str == "  ") return NO_PIECE;

    assert(std::tolower(str[0]) == 'a' || std::tolower(str[0]) == 'b');
    assert(std::tolower(str[1]) == 'o' || std::tolower(str[1]) == 'x');

    uint16_t piece = 0;

    if (std::tolower(str[0]) == 'b') setBit(piece, 0);
    if (std::tolower(str[1]) == 'x') setBit(piece, 1);
    if (std::isupper(str[0])) setBit(piece, 2);
    if (std::isupper(str[1])) setBit(piece, 3);

    return piece;
}

std::string cellToString(uint16_t cell)
{
    uint16_t row = cellToRow(cell);
    uint16_t col = cellToCol(cell);

    char first = 'a' + col;
    char second = '1' + (NUM_ROWS - row - 1);

    std::string str;
    str += first;
    str += second;

    return str;
}

uint16_t stringToCell(std::string str)
{
    assert(str.size() == 2);

    uint16_t row = NUM_ROWS - (str[1] - '1') - 1;
    uint16_t col = str[0] - 'a';

    assert(row >= 0 && row < NUM_ROWS);
    assert(col >= 0 && col < NUM_COLS);

    uint16_t cell = rowColToCell(row, col);

    return cell;
}

void play()
{
    State state;

    uint16_t player = 0;

    int currMove = 0;
    int minEvalMove = 15;

    while (true)
    {
        std::cout << "Board:" << std::endl;
        std::cout << "+----+----+----+----+";
        std::cout << std::endl;
        FOR_CELLS(i)
        {
            std::cout << "| ";
            std::cout << pieceToString(state.getPiece(i));
            std::cout << " ";

            uint16_t col = cellToCol(i);

            if (col == NUM_COLS - 1)
            {
                std::cout << "|";
                std::cout << std::endl;
                std::cout << "+----+----+----+----+";
                std::cout << std::endl;
            }
        }
        std::cout << std::endl;

        std::cout << "Player:" << std::endl;
        std::cout << player + 1 << std::endl;
        std::cout << std::endl;

        if (currMove++ >= minEvalMove)
        {
            std::cout << "Eval:" << std::endl;
            std::cout << eval(state) << std::endl;
            std::cout << std::endl;
        }

        if (state.isWon())
        {
            std::cout << "Win" << std::endl;
            break;
        }

        if (state.isDone())
        {
            std::cout << "Draw" << std::endl;
            break;
        }

        std::cout << "Piece:" << std::endl;
        std::string pieceStr;
        std::cin >> pieceStr;
        std::cout << pieceStr << std::endl;
        std::cout << std::endl;
        state.moveSelect(stringToPiece(pieceStr));

        player = 1 - player;

        std::cout << "Player:" << std::endl;
        std::cout << player + 1 << std::endl;
        std::cout << std::endl;

        if (currMove++ >= minEvalMove)
        {
            std::cout << "Eval:" << std::endl;
            std::cout << eval(state) << std::endl;
            std::cout << std::endl;
        }

        std::cout << "Cell:" << std::endl;
        std::string cellStr;
        std::cin >> cellStr;
        std::cout << cellStr << std::endl;
        std::cout << std::endl;
        state.movePlace(stringToCell(cellStr));
    }
}

void genRandGame(int seed)
{
    srand(seed);

    std::vector<std::string> pieces;
    FOR_PIECES(i)
    {
        pieces.push_back(pieceToString(i));
    }

    std::vector<std::string> cells;
    FOR_CELLS(i)
    {
        cells.push_back(cellToString(i));
    }

    std::random_shuffle(pieces.begin(), pieces.end());
    std::random_shuffle(cells.begin(), cells.end());

    for (uint16_t i = 0; i < std::min(NUM_PIECES, NUM_CELLS); ++i)
    {
        std::cout << pieces[i] << std::endl;
        std::cout << cells[i] << std::endl;
    }
}

int main()
{
    // genRandGame(42);

    play();

    return 0;
}
