
#include "fvCFD.H"
#include <vector>
#include <string>
#include <stdexcept>
#include <fstream>
struct fieldMeta
{
    word name;
    word className;
    fieldKind kind;
    label nCells;
    label nComp;
    fileName outFile;
};
enum class fieldKind
{
    SCALAR,
    VECTOR,
    SYMM_TENSOR,
    TENSOR,
    UNKNOWN
};
std::string shapeString(const std::vector<std::size_t>& shape)
{
    std::ostringstream os;
    os << "(";
    for (std::size_t i = 0; i < shape.size(); ++i)
    {
        os << shape[i];
        if (shape.size() == 1 || i + 1 < shape.size())
        {
            os << ",";
            if (i + 1 < shape.size())
            {
                os << " ";
            }
        }
    }
    os << ")";
    return os.str();
}

std::string makeNpyHeader
(
    const std::vector<std::size_t>& shape,
    const npyType dtype,
    const bool fortranOrder
)
{
    std::ostringstream dict;
    std::string npy_dtype = (dtype == npyType::FLOAT32 ? "<f4" : "<f8");
    dict
        << "{'descr': '" << npy_dtype
        << "', 'fortran_order': " << (fortranOrder ? "True" : "False")
        << ", 'shape': " << shapeString(shape)
        << ", }";

    std::string header = dict.str();

    const std::size_t preamble = 10;
    const std::size_t total = preamble + header.size() + 1;
    const std::size_t pad = (16 - (total % 16)) % 16;

    header.append(pad, ' ');
    header.push_back('\n');

    return header;
}

std::size_t writeNpyHeader
(
    std::ofstream& os,
    const std::vector<std::size_t>& shape,
    const npyType dtype,
    const bool fortranOrder
)
{
    const char magic[] = "\x93NUMPY";
    const unsigned char major = 1;
    const unsigned char minor = 0;

    const std::string header = makeNpyHeader(shape, dtype, fortranOrder);
    const std::uint16_t hlen = static_cast<std::uint16_t>(header.size());

    os.write(magic, 6);
    os.put(static_cast<char>(major));
    os.put(static_cast<char>(minor));

    char lenBytes[2];
    lenBytes[0] = static_cast<char>(hlen & 0xFF);
    lenBytes[1] = static_cast<char>((hlen >> 8) & 0xFF);
    os.write(lenBytes, 2);

    os.write(header.data(), header.size());

    return 10 + header.size();
}

Template <T>
class npyWriter
{
    std::ofstream os_;
    std::vector<std::size_t> shape_;
    npyType dtype_;
    bool fortranOrder_;
    metadata meta_;
    std::size_t dataStart_;
    const label nComp_;
    const label nTimes_;
    const label nCells_;
public:

    npyWriter
    (
        const metadata& meta,
        const std::vector<std::size_t>& shape,
        const npyType dtype,
        const bool fortranOrder

    )
    :
        meta_(meta),
        shape_(shape),
        dtype_(dtype),
        fortranOrder_(fortranOrder),
        dataStart_(0)
    {
        const fileName& fName = meta.outDir;
        nComp_ = meta.nComp;
        nCells_ = meta.nCells;
        nTimes_ = selectedTimes.size();

        os_.open(fName.c_str(), std::ios::binary | std::ios::trunc);

        if (!os_.good())
        {
            FatalErrorInFunction
                << "Cannot open output file " << fName
                << exit(FatalError);
        }

        dataStart_ = writeNpyHeader(os_, shape_, dtype_, fortranOrder_);
    }

    void seekToElement(std::size_t elem)
    {
        
        os_.seekp(dataStart_ + elem*dtypeSize(dtype_), std::ios::beg);
    }

    void write(const T fld, const timei)
    {


    }

    const label getElem(cell, comp, time)
    {
        if (fortranOrder_)
        {
            return 
        
        }
        else
        {
            return
        }
    }
    void flush()
    {
        os_.flush();
    }
};

