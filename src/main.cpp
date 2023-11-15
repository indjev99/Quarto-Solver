#include <iostream>
#include <cassert>
#include <string>
#include <cctype>
#include <array>
#include <vector>
#include <algorithm>
#include <map>
#include <fstream>
#include <random>
#include <ctime>

#define FOR_PROPS(i) for (uint16_t i = 0; i < NUM_PROPS; ++i)
#define FOR_PROPS_VARS(i, j) for (uint16_t i = 0; i < NUM_PROPS; ++i) for (uint16_t j = 0; j < NUM_VARS; ++j) 
#define FOR_PIECES(i) for (uint16_t i = 0; i < NUM_PIECES; ++i)
#define FOR_CELLS(i) for (uint16_t i = 0; i < NUM_CELLS; ++i)
#define FOR_WIN_LEN(i) for (uint16_t i = 0; i < WIN_LEN; ++i)

using int128_t = __int128;
using uint128_t = unsigned __int128;

std::mt19937 generator(time(nullptr));

template <class T>
auto randElem(const T& cont)
{
    std::uniform_int_distribution<uint32_t> distr(0, cont.size() - 1);
    return cont[distr(generator)];
}

constexpr uint16_t NUM_VARS = 2;
constexpr uint16_t NUM_PROPS = 4;
constexpr uint16_t NUM_PIECES = 1 << NUM_PROPS;
constexpr uint16_t NO_PIECE = NUM_PIECES;

constexpr uint16_t NUM_ROWS = 4;
constexpr uint16_t NUM_COLS = 4;
constexpr uint16_t NUM_CELLS = NUM_ROWS * NUM_COLS;

constexpr uint16_t NUM_MOVES = std::min(NUM_PIECES, NUM_CELLS);

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
    return (val >> n) & 1;
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

std::vector<uint16_t> computeRegWinMasks()
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

std::vector<uint16_t> computeTorusWinMasks()
{
    std::vector<uint16_t> winMasks;

    FOR_CELLS(i)
    {
        uint16_t row = cellToRow(i);
        uint16_t col = cellToCol(i);

        {
            uint16_t winMask = 0;
            FOR_WIN_LEN(j)
            {
                setBit(winMask, rowColToCell(row, (col + j) % NUM_COLS));
            }
            winMasks.push_back(winMask);
        }

        {
            uint16_t winMask = 0;
            FOR_WIN_LEN(j)
            {
                setBit(winMask, rowColToCell((row + j) % NUM_ROWS, col));
            }
            winMasks.push_back(winMask);
        }

        {
            uint16_t winMask = 0;
            FOR_WIN_LEN(j)
            {
                setBit(winMask, rowColToCell((row + j) % NUM_ROWS, (col + j) % NUM_COLS));
            }
            winMasks.push_back(winMask);
        }

        {
            uint16_t winMask = 0;
            FOR_WIN_LEN(j)
            {
                setBit(winMask, rowColToCell((row + j) % NUM_ROWS, (NUM_COLS + col - j) % NUM_COLS));
            }
            winMasks.push_back(winMask);
        }

        {
            uint16_t winMask = 0;
            FOR_WIN_LEN(j)
            {
                setBit(winMask, rowColToCell((row + j / WIN_SQ_SIDE) % NUM_ROWS, (col + j % WIN_SQ_SIDE) % NUM_COLS));
            }
            winMasks.push_back(winMask);
        }
    }

    return winMasks;
}

const std::vector<uint16_t> winMasks = computeTorusWinMasks();

std::array<std::vector<uint16_t>, NUM_CELLS> computeCellWinMasks()
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

const std::array<std::vector<uint16_t>, NUM_CELLS> cellWinMasks = computeCellWinMasks();

constexpr uint32_t NUM_CELL_MASKS = 1 << NUM_CELLS;

constexpr uint16_t NUM_ROTS = 8;

std::array<std::array<uint16_t, NUM_CELL_MASKS>, NUM_ROTS> computeRotMasks()
{
    static_assert(NUM_ROWS == NUM_COLS);

    std::array<std::array<uint16_t, NUM_CELL_MASKS>, NUM_ROTS> rotMasks;

    for (uint16_t rot = 0; rot < NUM_ROTS; ++rot)
    {
        for (uint32_t mask = 0; mask < NUM_CELL_MASKS; ++mask)
        {
            if (rot == 0)
            {
                rotMasks[rot][mask] = mask;
                continue;
            }

            uint16_t prev = rotMasks[rot != NUM_ROTS / 2 ? rot - 1 : 0][mask];
            uint16_t next = 0;

            FOR_CELLS(i)
            {
                uint16_t row = cellToRow(i);
                uint16_t col = cellToCol(i);
                uint16_t row2 = rot != NUM_ROTS / 2 ? NUM_COLS - col - 1 : NUM_COLS - row - 1;
                uint16_t col2 = rot != NUM_ROTS / 2 ? row : col;
                uint16_t i2 = rowColToCell(row2, col2);
                if (getBit(prev, i)) setBit(next, i2);
            }

            rotMasks[rot][mask] = next;
        }
    }

    return rotMasks;
}

const std::array<std::array<uint16_t, NUM_CELL_MASKS>, NUM_ROTS> rotMasks = computeRotMasks();

struct State
{
    uint16_t movesLeft;
    uint16_t currPiece;
    uint16_t piecesTaken;
    uint16_t cellsTaken;
    uint16_t cellsProps[NUM_PROPS][NUM_VARS];

    State()
    {
        movesLeft = NUM_MOVES;
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
        uint128_t minKey = -1;
    
        uint16_t minCellsTaken = -1;

        for (uint16_t rot = 0; rot < NUM_ROTS; ++rot)
        {
            minCellsTaken = std::min(minCellsTaken, rotMasks[rot][cellsTaken]);
        }

        for (uint16_t rot = 0; rot < NUM_ROTS; ++rot)
        {
            if (rotMasks[rot][cellsTaken] != minCellsTaken) continue;

            uint16_t otherCellsProps[NUM_PROPS];

            FOR_PROPS(i)
            {
                otherCellsProps[i] = std::min(rotMasks[rot][cellsProps[i][0]], rotMasks[rot][cellsProps[i][1]]);
            }

            std::sort(otherCellsProps, otherCellsProps + NUM_PROPS);

            uint128_t key = minCellsTaken;
            FOR_PROPS(i)
            {
                key = ((key << 16) | otherCellsProps[i]);
            }
    
            minKey = std::min<uint128_t>(minKey, key);
        }

        return minKey;
    }
};

constexpr uint16_t KEY_BITS = 56;
constexpr uint64_t KEY_BITS_MASK = (1ull << KEY_BITS) - 1;

struct TransTable
{
    struct Entry
    {
        uint64_t key : KEY_BITS = 0;
        int64_t val : 6 = 0;
        uint64_t isAlpha : 1 = 0;
        uint64_t isBeta : 1 = 0;
    };

    static_assert(sizeof(Entry) == 8);

    std::vector<Entry> data;

    TransTable(uint64_t size = 16782823): data(size) {}

    uint64_t index(uint128_t key) const
    {
        return key % data.size();
    }

    void clear()
    {
        std::fill(data.begin(), data.end(), Entry());
    }

    void put(uint128_t key, int16_t val, bool isAlpha, bool isBeta)
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
        if (data[i].key != (key & KEY_BITS_MASK)) return nullptr;
        return &data[i];
    }
};

constexpr uint16_t TT_DIV = 5;
constexpr uint16_t NUM_TTS = (NUM_MOVES - 2) / TT_DIV + 1;

TransTable transTables[NUM_TTS];

std::array<uint16_t, NUM_CELL_MASKS> computeLosePropVarCells()
{
    std::array<uint16_t, NUM_CELL_MASKS> losePropVarCells;

    uint16_t lastMask = NUM_CELL_MASKS - 1;

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

const std::array<uint16_t, NUM_CELL_MASKS> losePropVarCells = computeLosePropVarCells();

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

constexpr int16_t INF = (1 << 15) - 1;

uint64_t totalEvalStates;

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

            if ((state.cellsTaken & winMask) == winMask)
            {
                FOR_PROPS(j)
                {
                    prior += 2 * ((state.cellsProps[j][!getBit(piece, j)] & winMask) == winMask);
                }
            }
            else
            {
                FOR_PROPS(j)
                {
                    uint16_t leftover = (state.cellsProps[j][getBit(piece, j)] & winMask) ^ winMask;
                    prior -= (leftover & (leftover - 1)) == 0;
                }
            }
        }

        if (state.movesLeft == 2) prior = -prior;

        cellsPriors[numMoves++] = {prior, i};
    }

    std::sort(cellsPriors, cellsPriors + numMoves);

    int16_t val = -INF;

    for (auto ptr = cellsPriors; ptr != cellsPriors + numMoves; ++ptr)
    {
        uint16_t i = ptr->second;

        state.movePlace(i);
        int16_t nextVal = evalSelect(state, alpha, beta);
        state.undoPlace(piece, i);

        if (nextVal > val)
        {
            val = nextVal;
            alpha = std::max(alpha, val);
            if (alpha >= beta) break;
        }
    }

    assert(val > -INF);

    return val;
}

int16_t evalSelect(State& state, int16_t alpha, int16_t beta)
{
    ++totalEvalStates;

    assert(alpha < beta);

    assert(state.movesLeft > 0);
    assert(alpha >= - state.movesLeft);
    assert(beta <= state.movesLeft - 1);

    int16_t oldAlpha = alpha;
    int16_t oldBeta = beta;

    bool losePropsVars[NUM_PROPS][NUM_VARS];

    FOR_PROPS_VARS(i, j)
    {
        losePropsVars[i][j] = ~state.cellsTaken & losePropVarCells[state.cellsProps[i][j]];
    }

    uint16_t loseMask = 0;

    FOR_PROPS_VARS(i, j)
    {
        if (losePropsVars[i][j]) setBit(loseMask, i * NUM_VARS + j);
    }

    bool haveMove = false;

    for (uint16_t i : notLosingSelects[loseMask])
    {
        if (state.isPieceFree(i))
        {
            haveMove = true;
            break;
        }
    }

    if (!haveMove) return - state.movesLeft;

    if (state.movesLeft == 1) return 0;

    alpha = std::max<int16_t>(alpha, - (state.movesLeft - 2));

    if (alpha >= beta) return alpha;

    uint128_t key = state.getKey();

    uint16_t movesDone = NUM_MOVES - state.movesLeft;

    TransTable& currTransTable = transTables[movesDone / TT_DIV];

    const TransTable::Entry* entry = currTransTable.get(key);

    if (entry != nullptr && entry->isAlpha && entry->val > alpha)
    {
        alpha = entry->val;
        if (alpha >= beta) return alpha;
    }

    if (entry != nullptr && entry->isBeta && entry->val < beta)
    {
        beta = entry->val;
        if (alpha >= beta) return beta;
    }

    int16_t val = -INF;

    for (uint16_t i : notLosingSelects[loseMask])
    {
        if (!state.isPieceFree(i)) continue;

        state.moveSelect(i);
        int16_t nextVal = -evalPlace(state, -beta, -alpha);
        state.undoSelect();

        if (nextVal > val)
        {
            val = nextVal;
            alpha = std::max(alpha, val);
            if (alpha >= beta) break;
        }
    }

    assert(val > -INF);

    currTransTable.put(key, val, val > oldAlpha, val < oldBeta);

    return val;
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
    totalEvalStates = 0;

    if (state.isWon()) return state.movesLeft + 1;
    if (state.isDone()) return 0;
    if (state.isToPlace() && checkWinInOne(state)) return state.movesLeft;
    if (state.isToPlace() && state.movesLeft == 1) return 0;

    int16_t min = - (state.movesLeft - (state.isToPlace() ? 1 : 0));
    int16_t max = state.movesLeft - (state.isToPlace() ? 2 : 1);
    auto evalFunc = state.isToPlace() ? evalPlace : evalSelect;

    return evalFunc(state, min, max);
}

std::vector<std::pair<int16_t, std::vector<uint16_t>>> evalMoves(State state)
{
    if (state.isWon() || state.isDone()) return {};

    std::map<int16_t, std::vector<uint16_t>> valsMoves;

    if (state.isToPlace())
    {
        uint16_t piece = state.currPiece;
        FOR_CELLS(i)
        {
            if (!state.isCellFree(i)) continue;

            state.movePlace(i);
            int16_t val = eval(state);
            state.undoPlace(piece, i);

            valsMoves[val].push_back(i);
        }
    }
    else
    {
        FOR_PIECES(i)
        {
            if (!state.isPieceFree(i)) continue;

            state.moveSelect(i);
            int16_t val = - eval(state);
            state.undoSelect();

            valsMoves[val].push_back(i);
        }
    }

    std::vector<std::pair<int16_t, std::vector<uint16_t>>> valsMovesVec(valsMoves.begin(), valsMoves.end());
    std::reverse(valsMovesVec.begin(), valsMovesVec.end());

    return valsMovesVec;
}

void recInitTransTable(State& state, uint16_t depth)
{
    eval(state);

    if (depth == 0) return;
    if (state.isWon()) return;
    if (state.isDone()) return;

    if (state.isToPlace())
    {
        uint16_t piece = state.currPiece;

        FOR_CELLS(i)
        {
            if (!state.isCellFree(i)) continue;

            state.movePlace(i);
            recInitTransTable(state, depth - 1);
            state.undoPlace(piece, i);
        }
    }
    else
    {
        FOR_PIECES(i)
        {
            if (!state.isPieceFree(i)) continue;

            state.moveSelect(i);
            recInitTransTable(state, depth - 1);
            state.undoSelect();
        }
    }
}

void initTransTable()
{
    std::string fileName = "torus.tt";

    std::ifstream inFile(fileName, std::ios::in | std::ios::binary);

    if (inFile.good())
    {
        std::cerr << "Loading transposition table: " << fileName << std::endl;

        inFile.read((char*) transTables[0].data.data(), sizeof(TransTable::Entry) * transTables[0].data.size());

        return;
    }

    inFile.close();

    for (uint16_t depth = 0; depth <= 8; ++depth)
    {
        std::cerr << "Populating transposition table with depth: " << depth << std::endl;

        State state;
        recInitTransTable(state, depth);
    
        std::ofstream outFile(fileName, std::ios::out | std::ios::binary);
        outFile.write((char*) transTables[0].data.data(), sizeof(TransTable::Entry) * transTables[0].data.size());
    }
}

std::string evalToString(const State& state, int16_t val)
{
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
    char second = '1' + row;

    std::string str;
    str += first;
    str += second;

    return str;
}

uint16_t stringToCell(std::string str)
{
    assert(str.size() == 2);

    uint16_t row = str[1] - '1';
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

        std::cout << "Eval:" << std::endl;
        std::cout << evalToString(state, eval(state)) << std::endl;
        std::cout << std::endl;

        std::cout << "Moves:" << std::endl;
        for (const auto& [val, moves] : evalMoves(state))
        {
            for (uint16_t move : moves)
            {
                std::cout << pieceToString(move) << " ";
            }
            std::cout << ": " << evalToString(state, val) << std::endl;
        }
        std::cout << std::endl;

        std::cout << "Piece:" << std::endl;
        std::string pieceStr;
        // std::cin >> pieceStr;
        pieceStr = pieceToString(randElem(evalMoves(state)[0].second));
        std::cout << pieceStr << std::endl;
        std::cerr << pieceStr << std::endl;
        std::cout << std::endl;

        state.moveSelect(stringToPiece(pieceStr));

        player = 1 - player;

        std::cout << "Player:" << std::endl;
        std::cout << player + 1 << std::endl;
        std::cout << std::endl;

        std::cout << "Eval:" << std::endl;
        std::cout << evalToString(state, eval(state)) << std::endl;
        std::cout << std::endl;

        std::cout << "Moves:" << std::endl;
        for (const auto& [val, moves] : evalMoves(state))
        {
            for (uint16_t move : moves)
            {
                std::cout << cellToString(move) << " ";
            }
            std::cout << ": " << evalToString(state, val) << std::endl;
        }
        std::cout << std::endl;

        std::cout << "Cell:" << std::endl;
        std::string cellStr;
        // std::cin >> cellStr;
        cellStr = cellToString(randElem(evalMoves(state)[0].second));
        std::cout << cellStr << std::endl;
        std::cerr << cellStr << std::endl;
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
    initTransTable();

    play();

    return 0;
}
