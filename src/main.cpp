#include <iostream>
#include <cassert>
#include <string>
#include <cctype>

#define FOR_PROPS(i) for (uint16_t i = 0; i < NUM_PROPS; ++i)
#define FOR_PROPS_VARS(i, j) for (uint16_t i = 0; i < NUM_PROPS; ++i) for (uint16_t j = 0; j < NUM_VARS; ++j) 
#define FOR_PIECES(i) for (uint16_t i = 0; i < NUM_PIECES; ++i)
#define FOR_CELLS(i) for (uint16_t i = 0; i < NUM_CELLS; ++i)
#define FOR_ROWS_COLS(i, j) for (uint16_t i = 1; i <= NUM_ROWS; ++i) for (uint16_t j = 1; j <= NUM_COLS; ++j)

constexpr uint16_t NUM_VARS = 2u;
constexpr uint16_t NUM_PROPS = 4u;
constexpr uint16_t NUM_PIECES = 1 << NUM_PROPS;
constexpr uint16_t NO_PIECE = NUM_PIECES;

constexpr uint16_t NUM_ROWS = 4;
constexpr uint16_t NUM_COLS = 4;
constexpr uint16_t NUM_CELLS = NUM_ROWS * NUM_COLS;

static_assert(NUM_PIECES == NUM_CELLS);

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

struct State
{
    uint16_t currPiece;
    uint16_t piecesTaken;
    uint16_t cellsTaken;
    uint16_t cellsProps[NUM_PROPS][NUM_VARS];

    State()
    {
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
        assert(currPiece == NO_PIECE);

        setBit(piecesTaken, piece);
        currPiece = piece;
    }

    void undoSelect()
    {
        assert(currPiece != NO_PIECE);

        clearBit(piecesTaken, currPiece);
        currPiece = NO_PIECE;
    }

    void movePlace(uint16_t cell)
    {
        assert(currPiece != NO_PIECE);

        setBit(cellsTaken, cell);

        FOR_PROPS(i)
        {
            setBit(cellsProps[i][getBit(currPiece, i)], cell);
        }

        currPiece = NO_PIECE;
    }

    void undoPlace(uint16_t piece, uint16_t cell)
    {
        assert(currPiece == piece);

        clearBit(cellsTaken, cell);

        FOR_PROPS(i)
        {
            clearBit(cellsProps[i][getBit(piece, i)], cell);
        }

        currPiece = piece;
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

uint16_t rowColToCell(uint16_t row, uint16_t col)
{
    return (row - 1) * NUM_COLS + (col - 1);
}

uint16_t cellToRow(uint16_t cell)
{
    return cell / NUM_ROWS + 1;
}

uint16_t cellToCol(uint16_t cell)
{
    return cell % NUM_COLS + 1;
}

std::string pieceToString(uint16_t piece)
{
    if (piece == NO_PIECE) return "  ";

    char first = getBit(piece, 0) ? 'b' : 'a';
    char second = getBit(piece, 1) ? 'x' : 'o';
    if (getBit(piece, 2)) first = std::toupper(first);
    if (getBit(piece, 3)) second  = std::toupper(second );

    std::string pieceStr;
    pieceStr += first;
    pieceStr += second;
    return pieceStr;
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

void play()
{
    State state;

    while (true)
    {
        std::cout << "Board:" << std::endl;
        FOR_ROWS_COLS(i, j)
        {
            std::cout << pieceToString(state.getPiece(rowColToCell(i, j)));

            if (j < NUM_COLS) std::cout << " | ";
            else if (i < NUM_ROWS)
            {
                std::cout << std::endl;
                std::cout << "---+----+----+---";
                std::cout << std::endl;
            }
            else std::cout << std::endl;
        }
        std::cout << std::endl;

        if (state.currPiece == NO_PIECE)
        {
            std::cout << "Piece:" << std::endl;
            std::string pieceStr;
            std::cin >> pieceStr;
            std::cout << std::endl;
            state.moveSelect(stringToPiece(pieceStr));
            continue;
        }

        std::cout << "Piece:" << std::endl;
        std::cout << pieceToString(state.currPiece) << std::endl;
        std::cout << std::endl;

        std::cout << "Cell:" << std::endl;
        uint16_t row, col;
        std::cin >> row >> col;
        state.movePlace(rowColToCell(row, col));
        std::cout << std::endl;
    }
}

int main()
{
    play();

    return 0;
}
