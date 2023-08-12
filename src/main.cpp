#include <iostream>
#include <cassert>
#include <string>
#include <cctype>
#include <array>
#include <vector>
#include <algorithm>

#define FOR_PROPS(i) for (uint16_t i = 0; i < NUM_PROPS; ++i)
#define FOR_PROPS_VARS(i, j) for (uint16_t i = 0; i < NUM_PROPS; ++i) for (uint16_t j = 0; j < NUM_VARS; ++j) 
#define FOR_PIECES(i) for (uint16_t i = 0; i < NUM_PIECES; ++i)
#define FOR_CELLS(i) for (uint16_t i = 0; i < NUM_CELLS; ++i)
#define FOR_WIN_LEN(i) for (uint16_t i = 0; i < WIN_LEN; ++i)

using int128_t = __int128;
using uint128_t = unsigned __int128;

constexpr uint16_t NUM_VARS = 2;
constexpr uint16_t NUM_PROPS = 4;
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

std::array<std::vector<uint16_t>, NUM_CELLS> computeCellWinMasks(const std::vector<uint16_t>& winMasks)
{
    std::array<std::vector<uint16_t>, NUM_CELLS> cellWinMasks;

    FOR_CELLS(i)
    {
        for (uint16_t winMask : winMasks)
        {
            if (getBit(winMask, i)) cellWinMasks[i].push_back(winMask);
        }
    }

    return cellWinMasks;
}

const std::vector<uint16_t> winMasks = computeWinMasks();

const std::array<std::vector<uint16_t>, NUM_CELLS> cellWinMasks = computeCellWinMasks(winMasks);

struct State
{
    uint16_t movesLeft;
    uint16_t currPiece;
    uint16_t piecesTaken;
    uint16_t cellsTaken;
    uint16_t cellsProps[NUM_PROPS][NUM_VARS];

    State()
    {
        movesLeft = std::min(NUM_PIECES, NUM_CELLS);
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
    }

    void undoSelect()
    {
        assert(isToPlace());

        clearBit(piecesTaken, currPiece);
        currPiece = NO_PIECE;
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

    bool isToPlace() const
    {
        return currPiece != NO_PIECE;
    }

    bool isPieceFree(uint16_t piece) const
    {
        return !getBit(piecesTaken, piece);
    }

    bool isCellFree(uint16_t cell) const
    {
        return !getBit(cellsTaken, cell);
    }

    bool isWon() const
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

    bool isDone() const
    {
        return movesLeft == 0;
    }

    uint16_t getPiece(uint16_t cell) const
    {
        if (!getBit(cellsTaken, cell)) return NO_PIECE;
    
        uint16_t piece = 0;
        FOR_PROPS(i)
        {
            if (getBit(cellsProps[i][1], cell)) setBit(piece, i);
        }

        return piece;
    }

    uint128_t getKey() const
    {
        uint128_t key = 0;
        key = ((key << 16) | cellsTaken);
        FOR_PROPS(i)
        {
            key = ((key << 16) | cellsProps[i][0]);
        }
        return key;
    }
};

struct TransTable
{
    struct Entry
    {
        uint128_t key : 80 = 0;
        int128_t val : 16 = 0;
        uint128_t isAlpha : 1 = 0;
        uint128_t isBeta : 1 = 0;
    };

    static_assert(sizeof(Entry) == 16);

    std::vector<Entry> data;

    TransTable(uint64_t size = 8388593): data(size) {}

    uint64_t index(uint128_t key) const
    {
        return key % data.size();
    }

    void clear()
    {
        std::fill(data.begin(), data.end(), Entry());
    }

    void put(uint128_t key, uint16_t val, bool isAlpha, bool isBeta)
    {
        uint64_t i = index(key);
        data[i].key = key;
        data[i].val = val;
        data[i].isAlpha = isAlpha;
        data[i].isBeta = isBeta;
    }

    const Entry* get(uint128_t key) const
    {
        uint64_t i = index(key);
        if (data[i].key != key) return nullptr;
        return &data[i];
    }
};

uint64_t totalEvalStates;

TransTable transTable;

constexpr uint32_t NUM_PROP_VAR_MASKS = 1 << NUM_CELLS;

std::array<uint16_t, NUM_PROP_VAR_MASKS> computeLosePropVarCells()
{
    std::array<uint16_t, NUM_PROP_VAR_MASKS> losePropVarCells;

    uint16_t lastMask = NUM_PROP_VAR_MASKS - 1;

    for (uint16_t propVarMask = 0; ; ++propVarMask)
    {
        losePropVarCells[propVarMask] = 0;

        FOR_CELLS(i)
        {
            if (getBit(propVarMask, i)) continue;

            for (uint16_t winMask : cellWinMasks[i])
            {
                clearBit(winMask, i);

                if ((propVarMask & winMask) == winMask)
                {
                    setBit(losePropVarCells[propVarMask], i);
                    break;
                }
            }
        }

        if (propVarMask == lastMask) break;
    }

    return losePropVarCells;
}

const std::array<uint16_t, NUM_PROP_VAR_MASKS> losePropVarCells = computeLosePropVarCells();

constexpr uint16_t NUM_LOSE_MASKS = 1 << (NUM_PROPS * NUM_VARS);

std::array<std::vector<uint16_t>, NUM_LOSE_MASKS> computeNotLosingSelects()
{
    std::array<std::vector<uint16_t>, NUM_LOSE_MASKS> notLosingSelects;

    for (uint16_t loseMask = 0; loseMask < NUM_LOSE_MASKS; ++loseMask)
    {
        FOR_PIECES(i)
        {
            bool bad = false;

            FOR_PROPS_VARS(j, k)
            {
                if (getBit(loseMask, j * NUM_VARS + k) && getBit(i, j) == k)
                {
                    bad = true;
                    break;
                }
            }

            if (!bad) notLosingSelects[loseMask].push_back(i);
        }
    }

    return notLosingSelects;
}

const std::array<std::vector<uint16_t>, NUM_LOSE_MASKS> notLosingSelects = computeNotLosingSelects();

int16_t evalSelect(State& state, int16_t alpha, int16_t beta);

int16_t evalPlace(State& state, int16_t alpha, int16_t beta)
{
    ++totalEvalStates;

    assert(alpha < beta);

    assert(state.movesLeft > 0);
    assert(alpha >= - (state.movesLeft - 1));
    assert(beta <= std::max<int16_t>(state.movesLeft - 2, 0));

    uint16_t piece = state.currPiece;

    std::pair<int16_t, uint16_t> cellsPriors[NUM_CELLS];
    uint16_t numMoves = 0;

    FOR_CELLS(i)
    {
        if (!state.isCellFree(i)) continue;

        int16_t prior = 0;
        for (uint16_t winMask : cellWinMasks[i])
        {
            clearBit(winMask, i);
            FOR_PROPS(j)
            {
                prior += 16 * ((state.cellsProps[j][!getBit(piece, j)] & winMask) == winMask);
                uint16_t leftover = (state.cellsProps[j][getBit(piece, j)] & winMask) ^ winMask;
                prior -= (leftover & (leftover - 1)) == 0;
            }
        }

        cellsPriors[numMoves++] = {prior, i};
    }

    std::sort(cellsPriors, cellsPriors + numMoves);

    for (auto ptr = cellsPriors; ptr != cellsPriors + numMoves; ++ptr)
    {
        uint16_t i = ptr->second;

        state.movePlace(i);
        int16_t nextVal = evalSelect(state, alpha, beta);
        state.undoPlace(piece, i);

        assert(alpha <= nextVal);
        assert(nextVal <= beta);

        alpha = std::max(alpha, nextVal);

        if (alpha == beta) break;
    }

    return alpha;
}

int16_t evalSelect(State& state, int16_t alpha, int16_t beta)
{
    ++totalEvalStates;

    assert(alpha < beta);

    assert(state.movesLeft > 0);
    assert(alpha >= - state.movesLeft);
    assert(beta <= state.movesLeft - 1);

    int16_t oldAlpha = alpha;

    bool losePropsVars[NUM_PROPS][NUM_VARS];

    FOR_PROPS_VARS(i, j)
    {
        losePropsVars[i][j] = ~state.cellsTaken & losePropVarCells[state.cellsProps[i][j]];
    }

    FOR_PROPS(i)
    {
        if (std::all_of(losePropsVars[i], losePropsVars[i] + NUM_VARS, [](bool x){return x;}))
            return std::max<int16_t>(alpha, std::min<int16_t>(beta, -state.movesLeft));
    }

    alpha = std::max<int16_t>(alpha, std::min<int16_t>(beta, std::min<int16_t>(- (state.movesLeft - 2), 0)));

    assert(alpha <= beta);

    if (alpha == beta) return alpha;

    uint128_t key = state.getKey();

    const TransTable::Entry* entry = transTable.get(key);

    if (entry != nullptr && entry->isAlpha) alpha = std::max<int16_t>(alpha, std::min<int16_t>(beta, entry->val));
    if (entry != nullptr && entry->isBeta) beta = std::min<int16_t>(beta, std::max<int16_t>(alpha, entry->val));

    assert(alpha <= beta);

    if (alpha == beta) return alpha;

    uint16_t loseMask = 0;

    FOR_PROPS_VARS(i, j)
    {
        if (losePropsVars[i][j]) setBit(loseMask, i * NUM_VARS + j);
    }

    for (uint16_t i : notLosingSelects[loseMask])
    {
        if (!state.isPieceFree(i)) continue;

        state.moveSelect(i);
        int16_t nextVal = -evalPlace(state, -beta, -alpha);
        state.undoSelect();

        assert(alpha <= nextVal);
        assert(nextVal <= beta);

        alpha = std::max(alpha, nextVal);

        if (alpha == beta) break;
    }

    transTable.put(key, alpha, alpha > oldAlpha, alpha < beta);

    return alpha;
}

bool checkWinInOne(State& state)
{
    FOR_CELLS(i)
    {
        if (!state.isCellFree(i)) continue;

        for (uint16_t winMask : cellWinMasks[i])
        {
            clearBit(winMask, i);

            FOR_PROPS(j)
            {
                if ((state.cellsProps[j][getBit(state.currPiece, j)] & winMask) == winMask)
                    return true;
            }
        }
    }

    return false;
}

int16_t eval(State state)
{
    // transTable.clear();
    totalEvalStates = 0;

    if (state.isWon()) return state.movesLeft + 1;
    if (state.isDone()) return 0;
    if (state.isToPlace() && checkWinInOne(state)) return state.movesLeft;

    int16_t min = - (state.movesLeft - (state.isToPlace() ? 1 : 0));
    int16_t max = std::max<int16_t>(state.movesLeft - (state.isToPlace() ? 2 : 1), 0);
    auto evalFunc = state.isToPlace() ? evalPlace : evalSelect;

    return evalFunc(state, min, max);
}

std::string evalToString(const State& state, int16_t val)
{
    std::cerr << "nodes: " << totalEvalStates << std::endl;

    if (val == 0) return "Draw";

    std::string str = val > 0 ? "Win in " : "Loss in ";

    str += std::to_string(state.movesLeft + 1 - std::abs(val));

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
    int minEvalMove = 7;

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
            std::cout << evalToString(state, eval(state)) << std::endl;
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
            std::cout << evalToString(state, eval(state)) << std::endl;
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
