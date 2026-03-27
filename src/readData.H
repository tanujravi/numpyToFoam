
#include "fvCFD.H"
#include <vector>
#include <string>
#include <stdexcept>
#include <fstream>

struct NpyMeta
{
    fileName path;
    word file;
    word field;                   // e.g. "U"
    label procNo{-1};             // e.g. 0
    std::vector<size_t> shape;
    bool is_f8{false};
    bool is_f4{false};
    bool fortranOrder{false}; 
    std::streamoff dataOffset{0}; // where the raw array bytes start
};
static NpyMeta readNpyMeta(const fileName& filename)
{
    std::ifstream f(filename.c_str(), std::ios::binary);
    if (!f) throw std::runtime_error("Cannot open " + std::string(filename));

    char magic[6];
    f.read(magic, 6);

    char v[2];
    f.read(v, 2);

    unsigned short header_len;
    f.read(reinterpret_cast<char*>(&header_len), 2);

    std::string header(header_len, ' ');
    f.read(&header[0], header_len);
    NpyMeta meta;
    meta.path = filename;
    meta.file = filename.name();
    word base = filename.lessExt().name(); // e.g. "U_proc_0"

    label pos = base.rfind('_');          // last underscore
    meta.procNo = readLabel(base.substr(pos+1)); // "0"
    word left = base.substr(0, pos);      // "U_proc"
    label pos2 = left.rfind('_');         // underscore before proc
    meta.field = left.substr(0, pos2);    // "U"

    meta.is_f8 = header.find("<f8") != std::string::npos;
    meta.is_f4 = header.find("<f4") != std::string::npos;

    if (!meta.is_f8 && !meta.is_f4)
        throw std::runtime_error("Unsupported dtype in " + std::string(filename));
    // Parse Fortran order
    if (header.find("fortran_order") != std::string::npos)
    {
        if (header.find("True") != std::string::npos)
            meta.fortranOrder = true;
        else if (header.find("False") != std::string::npos)
            meta.fortranOrder = false;
        else
            throw std::runtime_error("Could not parse fortran_order in " + std::string(filename));
    }
    else
    {
        throw std::runtime_error("Missing fortran_order in header of " + std::string(filename));
    }

    auto p1 = header.find("(");
    auto p2 = header.find(")");
    std::string shape_str = header.substr(p1 + 1, p2 - p1 - 1);

    size_t start = 0;
    while (true)
    {
        auto comma = shape_str.find(",", start);
        std::string token = shape_str.substr(start, comma - start);

        // trim spaces
        token.erase(0, token.find_first_not_of(" \t"));
        token.erase(token.find_last_not_of(" \t") + 1);

        if (!token.empty())
            meta.shape.push_back(std::stoul(token));

        if (comma == std::string::npos) break;
        start = comma + 1;
    }

    meta.dataOffset = f.tellg(); // current position is start of raw data
    return meta;
}

static inline size_t prod(const std::vector<size_t>& v)
{
    size_t p = 1;
    for (auto x : v) p *= x;
    return p;
}

static inline void require(bool cond, const std::string& msg)
{
    if (!cond) throw std::runtime_error(msg);
}

static inline double readAtIndex(std::ifstream& f, const bool is_f8, const std::streamoff bytePos)
{
    f.seekg(bytePos, std::ios::beg);

    if (is_f8)
    {
        double v;
        f.read(reinterpret_cast<char*>(&v), sizeof(double));
        return v;
    }
    else
    {
        float v;
        f.read(reinterpret_cast<char*>(&v), sizeof(float));
        return static_cast<double>(v);
    }
}
static void readScalarSnapshot
(
    const NpyMeta& meta,
    const label timeIndex,
    scalarField& snapshot
)
{
    require(meta.shape.size() == 2, "Expected scalar npy with 2 dims");

    const label n0 = meta.shape[0];
    const label n1 = meta.shape[1];

    label nCells, nTimes;

    // Shape is still interpreted as (nCells, nTimes)
    nCells = n0;
    nTimes = n1;

    require(timeIndex >= 0 && timeIndex < nTimes, "timeIndex out of range");

    snapshot.setSize(nCells);

    std::ifstream f(meta.path.c_str(), std::ios::binary);
    require(bool(f), "Cannot open " + std::string(meta.path));

    const size_t elemSize = meta.is_f8 ? sizeof(double) : sizeof(float);

    for (label celli = 0; celli < nCells; ++celli)
    {
        std::streamoff idx;

        if (!meta.fortranOrder)
        {
            // C-order: (nCells, nTimes)
            // idx = cell*nTimes + time
            idx = static_cast<std::streamoff>(celli)*nTimes + timeIndex;
        }
        else
        {
            // Fortran-order: (nCells, nTimes)
            // idx = cell + time*nCells
            idx = static_cast<std::streamoff>(celli)
                + static_cast<std::streamoff>(timeIndex)*nCells;
        }

        snapshot[celli] = readAtIndex(f, meta.is_f8, meta.dataOffset + idx*elemSize);
    }
}


static void readVectorSnapshot
(
    const NpyMeta& meta,
    const label timeIndex,
    vectorField& snapshot
)
{
    require(meta.shape.size() == 3, "Expected vector npy with 3 dims");

    const label nCells = meta.shape[0];
    const label nComp  = meta.shape[1];
    const label nTimes = meta.shape[2];

    require(nComp == 3, "Vector middle dim must be 3");
    require(timeIndex >= 0 && timeIndex < nTimes, "timeIndex out of range");

    snapshot.setSize(nCells);

    std::ifstream f(meta.path.c_str(), std::ios::binary);
    require(bool(f), "Cannot open " + std::string(meta.path));

    const size_t elemSize = meta.is_f8 ? sizeof(double) : sizeof(float);

    for (label celli = 0; celli < nCells; ++celli)
    {
        double vx = 0.0, vy = 0.0, vz = 0.0;

        for (label comp = 0; comp < 3; ++comp)
        {
            std::streamoff idx;

            if (!meta.fortranOrder)
            {
                // C-order: (nCells, 3, nTimes)
                // idx = ((cell*3 + comp)*nTimes + time)
                idx =
                    (static_cast<std::streamoff>(celli)*3 + comp)*nTimes
                  + timeIndex;
            }
            else
            {
                // Fortran-order: (nCells, 3, nTimes)
                // idx = cell + nCells*(comp + 3*time)
                idx =
                    static_cast<std::streamoff>(celli)
                  + static_cast<std::streamoff>(nCells)*
                    (comp + static_cast<std::streamoff>(3)*timeIndex);
            }

            const double val = readAtIndex(f, meta.is_f8, meta.dataOffset + idx*elemSize);

            if (comp == 0) vx = val;
            else if (comp == 1) vy = val;
            else vz = val;
        }

        snapshot[celli] = vector(vx, vy, vz);
    }
}

static void readSymmTensorSnapshot
(
    const NpyMeta& meta,
    const label timeIndex,
    symmTensorField& snapshot
)
{
    require(meta.shape.size() == 3, "Expected symmetric tensor npy with 3 dims");

    const label nCells = meta.shape[0];
    const label nComp  = meta.shape[1];
    const label nTimes = meta.shape[2];

    require(nComp == 6, "SymmTensor middle dim must be 6");
    require(timeIndex >= 0 && timeIndex < nTimes, "timeIndex out of range");

    snapshot.setSize(nCells);

    std::ifstream f(meta.path.c_str(), std::ios::binary);
    require(bool(f), "Cannot open " + std::string(meta.path));

    const size_t elemSize = meta.is_f8 ? sizeof(double) : sizeof(float);

    for (label celli = 0; celli < nCells; ++celli)
    {
        double c[6];

        for (label comp = 0; comp < 6; ++comp)
        {
            std::streamoff idx;

            if (!meta.fortranOrder)
            {
                // C-order: (nCells, 6, nTimes)
                // idx = ((cell*6 + comp)*nTimes + time)
                idx =
                    (static_cast<std::streamoff>(celli)*6 + comp)*nTimes
                  + timeIndex;
            }
            else
            {
                // Fortran-order: (nCells, 6, nTimes)
                // idx = cell + nCells*(comp + 6*time)
                idx =
                    static_cast<std::streamoff>(celli)
                  + static_cast<std::streamoff>(nCells)*
                    (comp + static_cast<std::streamoff>(6)*timeIndex);
            }

            c[comp] = readAtIndex(f, meta.is_f8, meta.dataOffset + idx*elemSize);
        }

        snapshot[celli] = symmTensor
        (
            c[0], // xx
            c[1], // xy
            c[2], // xz
            c[3], // yy
            c[4], // yz
            c[5]  // zz
        );
    }
}

static void readTensorSnapshot
(
    const NpyMeta& meta,
    const label timeIndex,
    tensorField& snapshot
)
{
    require(meta.shape.size() == 3, "Expected tensor npy with 3 dims");

    const label nCells = meta.shape[0];
    const label nComp  = meta.shape[1];
    const label nTimes = meta.shape[2];

    require(nComp == 9, "Tensor middle dim must be 9");
    require(timeIndex >= 0 && timeIndex < nTimes, "timeIndex out of range");

    snapshot.setSize(nCells);

    std::ifstream f(meta.path.c_str(), std::ios::binary);
    require(bool(f), "Cannot open " + std::string(meta.path));

    const size_t elemSize = meta.is_f8 ? sizeof(double) : sizeof(float);

    for (label celli = 0; celli < nCells; ++celli)
    {
        double c[9];

        for (label comp = 0; comp < 9; ++comp)
        {
            std::streamoff idx;

            if (!meta.fortranOrder)
            {
                // C-order: (nCells, 9, nTimes)
                // idx = ((cell*9 + comp)*nTimes + time)
                idx =
                    (static_cast<std::streamoff>(celli)*9 + comp)*nTimes
                  + timeIndex;
            }
            else
            {
                // Fortran-order: (nCells, 9, nTimes)
                // idx = cell + nCells*(comp + 9*time)
                idx =
                    static_cast<std::streamoff>(celli)
                  + static_cast<std::streamoff>(nCells)*
                    (comp + static_cast<std::streamoff>(9)*timeIndex);
            }

            c[comp] = readAtIndex(f, meta.is_f8, meta.dataOffset + idx*elemSize);
        }

        snapshot[celli] = tensor
        (
            c[0], c[1], c[2],
            c[3], c[4], c[5],
            c[6], c[7], c[8]
        );
    }
}