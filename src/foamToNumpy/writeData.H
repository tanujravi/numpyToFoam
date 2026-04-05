#ifndef writeData_H
#define writeData_H

#include "fvCFD.H"
#include "Pstream.H"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdint>

namespace Foam
{

enum class fieldKind
{
    SCALAR,
    VECTOR,
    SYMM_TENSOR,
    TENSOR,
    UNKNOWN
};

enum class npyType
{
    FLOAT32,
    FLOAT64
};

struct fieldMeta
{
    word name;
    word className;
    fieldKind kind;
    label nCells;
    label nComp;
    fileName outFile;
};


//------------------------------------------------------------------------------
// small helpers
//------------------------------------------------------------------------------

inline npyType parseNpyType(const word& dataType)
{
    if (dataType == "float32")
    {
        return npyType::FLOAT32;
    }
    if (dataType == "float64")
    {
        return npyType::FLOAT64;
    }

    FatalErrorInFunction
        << "Unsupported dataType: " << dataType << nl
        << "Valid options are float32 and float64"
        << exit(FatalError);

    return npyType::FLOAT64;
}


inline bool parseFortranOrder(const word& storageOrder)
{
    if (storageOrder == "C" || storageOrder == "c")
    {
        return false;
    }
    if (storageOrder == "F" || storageOrder == "f")
    {
        return true;
    }

    FatalErrorInFunction
        << "Unsupported storageOrder: " << storageOrder << nl
        << "Valid options are C and F"
        << exit(FatalError);

    return false;
}


inline std::size_t dtypeSize(const npyType dtype)
{
    return (dtype == npyType::FLOAT32 ? sizeof(float) : sizeof(double));
}


inline std::string dtypeDescr(const npyType dtype)
{
    return (dtype == npyType::FLOAT32 ? "<f4" : "<f8");
}


inline std::string shapeString(const std::vector<std::size_t>& shape)
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


inline std::string makeNpyHeader
(
    const std::vector<std::size_t>& shape,
    const npyType dtype,
    const bool fortranOrder
)
{
    std::ostringstream dict;
    dict
        << "{'descr': '" << dtypeDescr(dtype)
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


inline std::size_t writeNpyHeader
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


inline std::vector<std::size_t> makeShape
(
    const label nTimes,
    const label nCells,
    const label nComp
)
{
    // shape = (nCells, nComp, nTimes)
    if (nComp == 1)
    {
        return
        {
            static_cast<std::size_t>(nCells),
            static_cast<std::size_t>(nTimes)
        };
    }

    return
    {
        static_cast<std::size_t>(nCells),
        static_cast<std::size_t>(nComp),
        static_cast<std::size_t>(nTimes)
    };
}


inline fieldKind classifyField(const word& cls)
{
    if (cls == volScalarField::typeName) return fieldKind::SCALAR;
    if (cls == volVectorField::typeName) return fieldKind::VECTOR;
    if (cls == volSymmTensorField::typeName) return fieldKind::SYMM_TENSOR;
    if (cls == volTensorField::typeName) return fieldKind::TENSOR;
    return fieldKind::UNKNOWN;
}


//------------------------------------------------------------------------------
// component extraction
//------------------------------------------------------------------------------

inline scalar componentValue(const scalar& v, const label comp)
{
    if (comp != 0)
    {
        FatalErrorInFunction
            << "Invalid scalar component index " << comp
            << exit(FatalError);
    }
    return v;
}


inline scalar componentValue(const vector& v, const label comp)
{
    return v[comp];
}


inline scalar componentValue(const symmTensor& v, const label comp)
{
    return v[comp];
}


inline scalar componentValue(const tensor& v, const label comp)
{
    return v[comp];
}


//------------------------------------------------------------------------------
// npyWriter
//------------------------------------------------------------------------------

template<class GeoField>
class npyWriter
{
    std::ofstream os_;
    std::vector<std::size_t> shape_;
    npyType dtype_;
    bool fortranOrder_;
    fieldMeta meta_;
    std::size_t dataStart_;
    label nComp_;
    label nTimes_;
    label nCells_;

    void writeNumber(const scalar val)
    {
        if (dtype_ == npyType::FLOAT32)
        {
            const float out = static_cast<float>(val);
            os_.write(reinterpret_cast<const char*>(&out), sizeof(float));
        }
        else
        {
            const double out = static_cast<double>(val);
            os_.write(reinterpret_cast<const char*>(&out), sizeof(double));
        }
    }

public:

    npyWriter
    (
        const fieldMeta& meta,
        const std::vector<std::size_t>& shape,
        const npyType dtype,
        const bool fortranOrder,
        const label nTimes
    )
    :
        shape_(shape),
        dtype_(dtype),
        fortranOrder_(fortranOrder),
        meta_(meta),
        dataStart_(0),
        nComp_(meta.nComp),
        nTimes_(nTimes),
        nCells_(meta.nCells)
    {
        os_.open(meta_.outFile.c_str(), std::ios::binary | std::ios::trunc);

        if (!os_.good())
        {
            FatalErrorInFunction
                << "Cannot open output file " << meta_.outFile
                << exit(FatalError);
        }

        dataStart_ = writeNpyHeader(os_, shape_, dtype_, fortranOrder_);
    }


    void seekToElement(const std::size_t elem)
    {
        os_.seekp(dataStart_ + elem*dtypeSize(dtype_), std::ios::beg);
    }


    std::size_t getElem
    (
        const label celli,
        const label compi,
        const label timei
    ) const
    {
        if (nComp_ == 1)
        {
            // shape = (nCells, nTimes)
            if (fortranOrder_)
            {
                // Fortran-order: offset = cell + nCells*time
                return
                    static_cast<std::size_t>(celli)
                  + static_cast<std::size_t>(nCells_)*timei;
            }
            else
            {
                // C-order: offset = cell*nTimes + time
                return
                    static_cast<std::size_t>(celli)*nTimes_
                  + static_cast<std::size_t>(timei);
            }
        }
        else
        {
            // shape = (nCells, nComp, nTimes)

            if (fortranOrder_)
            {
                // offset = cell + nCells*(comp + nComp*time)
                return
                    static_cast<std::size_t>(celli)
                  + static_cast<std::size_t>(nCells_)*
                    (
                        static_cast<std::size_t>(compi)
                      + static_cast<std::size_t>(nComp_)*timei
                    );
            }
            else
            {
                // offset = ((cell*nComp) + comp)*nTimes + time
                return
                    (
                        static_cast<std::size_t>(celli)*nComp_
                      + static_cast<std::size_t>(compi)
                    )*nTimes_
                  + static_cast<std::size_t>(timei);
            }
        }
    }


    void write(const GeoField& fld, const label timei)
    {
        if (fld.size() != nCells_)
        {
            FatalErrorInFunction
                << "Field size mismatch for " << meta_.name
                << " at time index " << timei << nl
                << "Expected " << nCells_ << " cells but got " << fld.size()
                << exit(FatalError);
        }

        if (nComp_ == 1)
        {
            for (label celli = 0; celli < nCells_; ++celli)
            {
                const std::size_t elem = getElem(celli, 0, timei);
                seekToElement(elem);
                writeNumber(componentValue(fld[celli], 0));
            }
        }
        else
        {
            for (label celli = 0; celli < nCells_; ++celli)
            {
                for (label compi = 0; compi < nComp_; ++compi)
                {
                    const std::size_t elem = getElem(celli, compi, timei);
                    seekToElement(elem);
                    writeNumber(componentValue(fld[celli], compi));
                }
            }
        }
    }


    void flush()
    {
        os_.flush();
    }
};

} // End namespace Foam

#endif