
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
    require(meta.shape.size() == 2, "Expected scalar npy with 2 dims (nCells, nTimes)");
    const label nCells = meta.shape[0];
    const label nTimes = meta.shape[1];
    require(timeIndex >= 0 && timeIndex < nTimes, "timeIndex out of range");

    snapshot.setSize(nCells);

    std::ifstream f(meta.path.c_str(), std::ios::binary);
    require(bool(f), "Cannot open " + std::string(meta.path));

    const size_t elemSize = meta.is_f8 ? sizeof(double) : sizeof(float);

    for (label celli = 0; celli < nCells; ++celli)
    {
        const std::streamoff idx = static_cast<std::streamoff>(celli)*nTimes + timeIndex;
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
    require(meta.shape.size() == 3, "Expected vector npy with 3 dims (nCells, 3, nTimes)");
    const label nCells = meta.shape[0];
    const label nComp  = meta.shape[1];
    const label nTimes = meta.shape[2];
    require(nComp == 3, "Vector middle dim must be 3");
    require(timeIndex >= 0 && timeIndex < nTimes, "timeIndex out of range");

    snapshot.setSize(nCells);

    std::ifstream f(meta.path.c_str(), std::ios::binary);
    require(bool(f), "Cannot open " + std::string(meta.path));

    const size_t elemSize = meta.is_f8 ? sizeof(double) : sizeof(float);

    // row-major: idx = ((cell * 3 + comp) * nTimes + time)
    for (label celli = 0; celli < nCells; ++celli)
    {
        double vx = 0, vy = 0, vz = 0;
        {
            const std::streamoff idx = (static_cast<std::streamoff>(celli)*3 + 0)*nTimes + timeIndex;
            vx = readAtIndex(f, meta.is_f8, meta.dataOffset + idx*elemSize);
        }
        {
            const std::streamoff idx = (static_cast<std::streamoff>(celli)*3 + 1)*nTimes + timeIndex;
            vy = readAtIndex(f, meta.is_f8, meta.dataOffset + idx*elemSize);
        }
        {
            const std::streamoff idx = (static_cast<std::streamoff>(celli)*3 + 2)*nTimes + timeIndex;
            vz = readAtIndex(f, meta.is_f8, meta.dataOffset + idx*elemSize);
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
    require(meta.shape.size() == 4, "Expected tensor npy with 4 dims (nCells, 3, 3, nTimes)");
    const label nCells = meta.shape[0];
    const label nA     = meta.shape[1];
    const label nB     = meta.shape[2];
    const label nTimes = meta.shape[3];

    require(nA == 3 && nB == 3, "Tensor dims must be 3x3");
    require(timeIndex >= 0 && timeIndex < nTimes, "timeIndex out of range");

    snapshot.setSize(nCells);

    std::ifstream f(meta.path.c_str(), std::ios::binary);
    require(bool(f), "Cannot open " + std::string(meta.path));

    const size_t elemSize = meta.is_f8 ? sizeof(double) : sizeof(float);

    // row-major: idx = (((cell*3 + i)*3 + j)*nTimes + time)
    for (label celli = 0; celli < nCells; ++celli)
    {
        double a[3][3];

        for (label i = 0; i < 3; ++i)
        for (label j = 0; j < 3; ++j)
        {
            const std::streamoff idx =
                (((static_cast<std::streamoff>(celli)*3 + i)*3 + j)*nTimes + timeIndex);

            a[i][j] = readAtIndex(f, meta.is_f8, meta.dataOffset + idx*elemSize);
        }

        // Enforce symmetry (safe if stored tensor is nearly symmetric)
        const scalar xx = a[0][0];
        const scalar yy = a[1][1];
        const scalar zz = a[2][2];

        const scalar xy = 0.5*(a[0][1] + a[1][0]);
        const scalar xz = 0.5*(a[0][2] + a[2][0]);
        const scalar yz = 0.5*(a[1][2] + a[2][1]);

        // OpenFOAM symmTensor ctor order: xx, xy, xz, yy, yz, zz
        snapshot[celli] = symmTensor(xx, xy, xz, yy, yz, zz);
    }
}
