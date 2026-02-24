/*
 * Copyright 2018 Swinburne University of Technology
 * This file is part of TAO-ScienceModules.
 */

#ifndef anyBuf_h
#define anyBuf_h

#include <stdint.h>

#if !defined(__BIG_ENDIAN__) && (defined(WIN32) || defined(MIPSEL) || defined(LINUX) || \
                                 defined(LIN) || defined(i386) || defined(__x86_64))
#define bigEndian() (false)
#else
#define bigEndian() (true)
#endif // __BIG_ENDIAN__

class anyBuf
{
public:
    anyBuf()
    {
        int size = 1;

        buf = new uint8_t[size + 1];
        max_size = size;
        cur_size = 0;
        cur_offset = 0;
    }

    anyBuf(int32_t size)
    {
        buf = new uint8_t[size + 1];
        max_size = size;
        cur_size = 0;
        cur_offset = 0;
    }

    anyBuf(anyBuf& from)
    {
        max_size = from.max_size;
        cur_size = from.cur_size;
        cur_offset = from.cur_offset;
        buf = new uint8_t[from.max_size + 1];

        if (max_size > 0)
        {
            bcopy((void*)from.buf, (void*)buf, max_size);
        }
    }

    anyBuf& operator=(const anyBuf& from)
    {
        max_size = from.max_size;
        cur_size = from.cur_size;
        cur_offset = from.cur_offset;
        buf = new uint8_t[from.max_size + 1];

        if (max_size > 0)
        {
            bcopy((void*)from.buf, (void*)buf, max_size);
        }

        return *this;
    }

    virtual ~anyBuf()
    {
        delete[] buf;
        buf = NULL;
    }

    void clear()
    {
        cur_offset = 0;
        memset(buf, 0, max_size);
    }

    void seek(int32_t offset) { cur_offset = offset; }

    int32_t tell() { return cur_offset; }

    int32_t resize(int32_t size)
    {
        int32_t status = true;

        if (size > max_size)
        {
            uint8_t* newbuf = new uint8_t[size];

            if (newbuf == NULL)
            {
                return false;
            }
            int32_t max = cur_offset;
            if (cur_size > max)
                max = cur_size;

            if (max > 0)
            {
                bcopy((void*)buf, (void*)newbuf, max);
            }

            delete[] buf;
            buf = newbuf;
            max_size = size;
        }

        cur_size = size;
        return status;
    }

    int32_t getMaxSize() { return max_size; }

    int32_t getCurSize() { return cur_size; }

    int32_t getOffset() { return cur_offset; }

    uint8_t* getbuf() { return buf; }

    void skip(int32_t nskip) { cur_offset += nskip; }

    int32_t toBIL(int32_t bands, int32_t bytesPerElement)
    {
        int32_t i, j, k;
        uint8_t* oldbuf = new uint8_t[cur_offset];
        int32_t bytesPerBand = cur_offset / bands;
        int32_t elements = bytesPerBand / bytesPerElement;

        for (i = 0; i < cur_offset; i++)
        {
            oldbuf[i] = buf[i];
        }

        int32_t to = 0;
        int32_t from = 0;

        for (i = 0; i < elements; i++)
        {
            for (k = 0; k < bands; k++)
            {
                from = i * bands * bytesPerElement + k * bytesPerElement;
                to = k * bytesPerBand + i * bytesPerElement;

                for (j = 0; j < bytesPerElement; j++)
                {
                    buf[to++] = oldbuf[from++];
                }
            }
        }

        delete[] oldbuf;
        return 0;
    }

    int32_t fromBIL(int32_t bands, int32_t bytesPerElement)
    {
        int32_t i, j, k;
        uint8_t* oldbuf = new uint8_t[cur_offset];
        int32_t bytesPerBand = cur_offset / bands;
        int32_t elements = bytesPerBand / bytesPerElement;

        for (i = 0; i < cur_offset; i++)
        {
            oldbuf[i] = buf[i];
        }

        int32_t to = 0;
        int32_t from = 0;

        for (i = 0; i < elements; i++)
        {
            for (k = 0; k < bands; k++)
            {
                from = i * bands * bytesPerElement + k * bytesPerElement;
                to = k * bytesPerBand + i * bytesPerElement;

                for (j = 0; j < bytesPerElement; j++)
                {
                    buf[from++] = oldbuf[to++];
                }
            }
        }

        delete[] oldbuf;
        return 0;
    }

    uint8_t popUInt8()
    {

        uint8_t v = buf[cur_offset];
        cur_offset += 1;
        return v;
    }

    int8_t popInt8()
    {
        int8_t v;
        pop(1, (void*)&v, !bigEndian());
        return v;
    }

    int16_t popBigEndianInt16()
    {
        int16_t v;
        pop(2, (void*)&v, !bigEndian());
        return v;
    }

    uint16_t popBigEndianUInt16()
    {
        uint16_t v;
        pop(2, (void*)&v, !bigEndian());
        return v;
    }

    int32_t popBigEndianInt32()
    {
        int32_t v;
        pop(4, (void*)&v, !bigEndian());
        return v;
    }

    int64_t popBigEndianInt64()
    {
        int64_t v;
        pop(8, (void*)&v, !bigEndian());
        return v;
    }

    uint64_t popBigEndianU64()
    {
        uint64_t v;
        pop(8, (void*)&v, !bigEndian());
        return v;
    }

    uint32_t popBigEndianUInt32()
    {
        uint32_t v;
        pop(4, (void*)&v, !bigEndian());
        return v;
    }

    float popBigEndianFlt32()
    {
        float v;
        pop(4, (void*)&v, !bigEndian());
        return v;
    }

    double popBigEndianFlt64()
    {
        double v;
        pop(8, (void*)&v, !bigEndian());
        return v;
    }

    int32_t popByteArray(int32_t elements, uint8_t* values, int32_t skip_elements = 0)
    {
        int32_t bytesize = sizeof(uint8_t);
        for (int32_t i = 0; i < elements; i++)
        {
            pop(bytesize, &values[i], false);

            if (skip_elements > 0)
            {
                skip(skip_elements * bytesize);
            }
        }
        return 0;
    }

    int32_t popS08Array(int32_t elements, int8_t* values, int32_t skip_elements = 0)
    {
        int32_t bytesize = sizeof(int8_t);
        for (int32_t i = 0; i < elements; i++)
        {
            pop(bytesize, &values[i], false);

            if (skip_elements > 0)
            {
                skip(skip_elements * bytesize);
            }
        }
        return 0;
    }

    int32_t popBigEndianInt16Array(int32_t elements, int16_t* values, int32_t skip_elements = 0)
    {
        int32_t bytesize = sizeof(int16_t);
        for (int32_t i = 0; i < elements; i++)
        {
            pop(bytesize, &values[i], !bigEndian());

            if (skip_elements > 0)
            {
                skip(skip_elements * bytesize);
            }
        }
        return 0;
    }

    int32_t popBigEndianUInt16Array(int32_t elements, uint16_t* values, int32_t skip_elements = 0)
    {
        int32_t bytesize = sizeof(uint16_t);
        for (int32_t i = 0; i < elements; i++)
        {
            pop(bytesize, &values[i], !bigEndian());

            if (skip_elements > 0)
            {
                skip(skip_elements * bytesize);
            }
        }
        return 0;
    }

    int32_t popBigEndianInt32Array(int32_t elements, int32_t* values, int32_t skip_elements = 0)
    {
        int32_t bytesize = sizeof(int32_t);
        for (int32_t i = 0; i < elements; i++)
        {
            pop(bytesize, &values[i], !bigEndian());

            if (skip_elements > 0)
            {
                skip(skip_elements * bytesize);
            }
        }
        return 0;
    }

    int32_t popBigEndianUInt32Array(int32_t elements, uint32_t* values, int32_t skip_elements = 0)
    {
        int32_t bytesize = sizeof(uint32_t);
        for (int32_t i = 0; i < elements; i++)
        {
            pop(bytesize, &values[i], !bigEndian());

            if (skip_elements > 0)
            {
                skip(skip_elements * bytesize);
            }
        }
        return 0;
    }

    int32_t popBigEndianS64Array(int32_t elements, int64_t* values, int32_t skip_elements = 0)
    {
        int32_t bytesize = sizeof(int64_t);
        for (int32_t i = 0; i < elements; i++)
        {
            pop(bytesize, &values[i], !bigEndian());

            if (skip_elements > 0)
            {
                skip(skip_elements * bytesize);
            }
        }
        return 0;
    }

    int32_t popBigEndianU64Array(int32_t elements, uint64_t* values, int32_t skip_elements = 0)
    {
        int32_t bytesize = sizeof(uint64_t);
        for (int32_t i = 0; i < elements; i++)
        {
            pop(bytesize, &values[i], !bigEndian());

            if (skip_elements > 0)
            {
                skip(skip_elements * bytesize);
            }
        }
        return 0;
    }

    int32_t popBigEndianFlt32Array(int32_t elements, float* values, int32_t skip_elements = 0)
    {
        int32_t bytesize = sizeof(float);
        for (int32_t i = 0; i < elements; i++)
        {
            pop(bytesize, &values[i], !bigEndian());

            if (skip_elements > 0)
            {
                skip(skip_elements * bytesize);
            }
        }
        return 0;
    }

    int32_t popBigEndianFlt64Array(int32_t elements, double* values, int32_t skip_elements = 0)
    {
        int32_t bytesize = sizeof(double);
        for (int32_t i = 0; i < elements; i++)
        {
            pop(bytesize, &values[i], !bigEndian());

            if (skip_elements > 0)
            {
                skip(skip_elements * bytesize);
            }
        }
        return 0;
    }

    int16_t popLittleEndianInt16()
    {
        int16_t v;
        pop(sizeof(int16_t), (void*)&v, bigEndian());
        return v;
    }

    uint16_t popLittleEndianUInt16()
    {
        uint16_t v;
        pop(sizeof(uint16_t), (void*)&v, bigEndian());
        return v;
    }

    int32_t popLittleEndianInt32()
    {
        int32_t v;
        pop(sizeof(int32_t), (void*)&v, bigEndian());
        return v;
    }

    uint32_t popLittleEndianUInt32()
    {
        uint32_t v;
        pop(sizeof(uint32_t), (void*)&v, bigEndian());
        return v;
    }

    int64_t popLittleEndianS64()
    {
        int64_t v;
        pop(sizeof(int64_t), (void*)&v, bigEndian());
        return v;
    }

    uint64_t popLittleEndianU64()
    {
        uint64_t v;
        pop(sizeof(uint64_t), (void*)&v, bigEndian());
        return v;
    }

    float popLittleEndianFlt32()
    {
        float v;
        pop(sizeof(float), (void*)&v, bigEndian());
        return v;
    }

    double popLittleEndianFlt64()
    {
        double v;
        pop(sizeof(double), (void*)&v, bigEndian());
        return v;
    }

    int32_t popLittleEndianInt16Array(int32_t elements, int16_t* values, int32_t skip_elements = 0)
    {

        int32_t bytesize = sizeof(int16_t);
        for (int32_t i = 0; i < elements; i++)
        {
            pop(bytesize, &values[i], bigEndian());

            if (skip_elements > 0)
            {
                skip(skip_elements * bytesize);
            }
        }

        return 0;
    }

    int32_t popLittleEndianUInt16Array(int32_t elements, uint16_t* values,
                                       int32_t skip_elements = 0)
    {

        int32_t bytesize = sizeof(uint16_t);
        for (int32_t i = 0; i < elements; i++)
        {
            pop(bytesize, &values[i], bigEndian());

            if (skip_elements > 0)
            {
                skip(skip_elements * bytesize);
            }
        }

        return 0;
    }

    int32_t popLittleEndianInt32Array(int32_t elements, int32_t* values, int32_t skip_elements = 0)
    {

        int32_t bytesize = sizeof(int32_t);
        for (int32_t i = 0; i < elements; i++)
        {
            pop(bytesize, &values[i], bigEndian());

            if (skip_elements > 0)
            {
                skip(skip_elements * bytesize);
            }
        }

        return 0;
    }

    int32_t popLittleEndianUInt32Array(int32_t elements, uint32_t* values,
                                       int32_t skip_elements = 0)
    {

        int32_t bytesize = sizeof(uint32_t);
        for (int32_t i = 0; i < elements; i++)
        {
            pop(bytesize, &values[i], bigEndian());

            if (skip_elements > 0)
            {
                skip(skip_elements * bytesize);
            }
        }

        return 0;
    }

    int32_t popLittleEndianS64Array(int32_t elements, int64_t* values, int32_t skip_elements = 0)
    {

        int32_t bytesize = sizeof(int64_t);
        for (int32_t i = 0; i < elements; i++)
        {
            pop(bytesize, &values[i], bigEndian());

            if (skip_elements > 0)
            {
                skip(skip_elements * bytesize);
            }
        }

        return 0;
    }

    int32_t popLittleEndianU64Array(int32_t elements, uint64_t* values, int32_t skip_elements = 0)
    {

        int32_t bytesize = sizeof(uint64_t);
        for (int32_t i = 0; i < elements; i++)
        {
            pop(bytesize, &values[i], bigEndian());

            if (skip_elements > 0)
            {
                skip(skip_elements * bytesize);
            }
        }

        return 0;
    }

    int32_t popLittleEndianFlt32Array(int32_t elements, float* values, int32_t skip_elements = 0)
    {

        int32_t bytesize = sizeof(float);
        for (int32_t i = 0; i < elements; i++)
        {
            pop(bytesize, &values[i], bigEndian());

            if (skip_elements > 0)
            {
                skip(skip_elements * bytesize);
            }
        }

        return 0;
    }

    int32_t popLittleEndianFlt64Array(int32_t elements, double* values, int32_t skip_elements = 0)
    {

        int32_t bytesize = sizeof(double);
        for (int32_t i = 0; i < elements; i++)
        {
            pop(bytesize, &values[i], bigEndian());

            if (skip_elements > 0)
            {
                skip(skip_elements * bytesize);
            }
        }

        return 0;
    }
    /*

    Int32 pushByte(const Byte v);

    Int32 pushS08(const s08 v);

    Int32 pushBigEndianInt16(const Int16 v);

    Int32 pushBigEndianUInt16(const UInt16 v);

    Int32 pushBigEndianInt32(const Int32 v);

    Int32 pushBigEndianLogical(const Logical v);

    Int32 pushBigEndianUInt32(const UInt32 v);

    Int32 pushBigEndianS64(const s64 v);

    Int32 pushBigEndianU64(const u64 v);

    Int32 pushBigEndianFlt32(const Flt32 v);

    Int32 pushBigEndianFlt64(const Flt64 v);

    Int32 pushBigEndianVector3lf(const Vector3lf& v);

    Int32 pushBigEndianTensor(const Tensor& v);

    Int32 pushLittleEndianLogical(const Logical v);

    Int32 pushLittleEndianInt16(const Int16 v);

    Int32 pushLittleEndianUInt16(const UInt16 v);

    Int32 pushLittleEndianInt24(const Int32 v);

    Int32 pushLittleEndianInt32(const Int32 v);

    Int32 pushLittleEndianS64(const s64 v);

    Int32 pushLittleEndianU64(const u64 v);

    Int32 pushLittleEndianFlt32(const Flt32 v);

    Int32 pushLittleEndianFlt64(const Flt64 v);

    Int32 pushLittleEndianVector3lf(const Vector3lf& v);

    Int32 pushLittleEndianTensor(const Tensor& v);

    Int32 pushIbmFlt32(const Flt32 v);

    Int32 pushCString(Byte *inbuf, Int32 maxlength);

    Int32 pushCString(dfastr inbuf, Int32 maxlength);

    Int32 pushVString(dfastr inbuf);

    static Flt32 ibm2ieee(const Flt32 ibm);

    static Flt32 ieee2ibm(const Flt32 ieee);

    Int64 getBits(Int32 frombit, Int32 nbits);
     */

    inline void swap4(char* data, int num)
    {
        long k, i, j, l;
        char tmp;

        for (k = 0; k < num; k++)
        {
            l = k * 4;

            for (i = l, j = l + 3; i < j; i++, j--)
            {
                tmp = data[i];
                data[i] = data[j];
                data[j] = tmp;
            }
        }
    }

    /*
    Int32 anyBuf::pushByte(const Byte iv)
    {
        assert(cur_offset < max_size && "push past bounds");
        buf[cur_offset] = iv;
        cur_offset     += 1;
        return cur_offset;
    }


    Int32 anyBuf::pushS08(const s08 iv)
    {
        bcopy((void*) &iv, (void*) &buf[cur_offset], 1);
        cur_offset += 1;
        return cur_offset;
    }


    Int32 anyBuf::pushBigEndianInt16(const Int16 iv)
    {
        Int16 v = iv;

        if (!bigEndian())
        {
            swap_bytes((char*) &v, 2);
        }

        bcopy((void*) &v, (void*) &buf[cur_offset], 2);
        cur_offset += 2;
        return cur_offset;
    }


    Int32 anyBuf::pushBigEndianUInt16(const UInt16 iv)
    {
        UInt16 v = iv;

        if (!bigEndian())
        {
            swap_bytes((char*) &v, 2);
        }

        bcopy((void*) &v, (void*) &buf[cur_offset], 2);
        cur_offset += 2;
        return cur_offset;
    }


    Int32 anyBuf::pushBigEndianInt32(const Int32 iv)
    {
        Int32 v = iv;

        if (!bigEndian())
        {
            swap_bytes((char*) &v, 4);
        }

        bcopy((void*) &v, (void*) &buf[cur_offset], 4);
        cur_offset += 4;
        return cur_offset;
    }

    Int32 anyBuf::pushBigEndianLogical(const Logical iv)
    {
        Logical v = iv;

        if (!bigEndian())
        {
            swap_bytes((char*) &v, sizeof(Logical));
        }

        bcopy((void*) &v, (void*) &buf[cur_offset], sizeof(Logical));
        cur_offset += sizeof(Logical);
        return cur_offset;
    }


    Int32 anyBuf::pushBigEndianUInt32(const UInt32 iv)
    {
        UInt32 v = iv;

        if (!bigEndian())
        {
            swap_bytes((char*) &v, 4);
        }

        bcopy((void*) &v, (void*) &buf[cur_offset], 4);
        cur_offset += 4;
        return cur_offset;
    }


    Int32 anyBuf::pushBigEndianS64(const s64 iv)
    {
        s64 v = iv;

        if (!bigEndian())
        {
            swap_bytes((char*) &v, 8);
        }

        bcopy((void*) &v, (void*) &buf[cur_offset], 8);
        cur_offset += 8;
        return cur_offset;
    }


    Int32 anyBuf::pushBigEndianU64(const u64 iv)
    {
        u64 v = iv;

        if (!bigEndian())
        {
            swap_bytes((char*) &v, 8);
        }

        bcopy((void*) &v, (void*) &buf[cur_offset], 8);
        cur_offset += 8;
        return cur_offset;
    }


    Int32 anyBuf::pushBigEndianFlt32(const Flt32 iv)
    {
        Flt32 v = iv;

        if (!bigEndian())
        {
            swap_bytes((char*) &v, 4);
        }

        bcopy((void*) &v, (void*) &buf[cur_offset], 4);
        cur_offset += 4;
        return cur_offset;
    }


    Int32 anyBuf::pushBigEndianFlt64(const Flt64 iv)
    {
        Flt64 v = iv;

        if (!bigEndian())
        {
            swap_bytes((char*) &v, 8);
        }

        bcopy((void*) &v, (void*) &buf[cur_offset], 8);
        cur_offset += 8;
        return cur_offset;
    }

    Int32 anyBuf::pushBigEndianVector3lf(const Vector3lf& iv)
    {
        double x = iv.X();
        double y = iv.Y();
        double z = iv.Z();

        if (!bigEndian())
        {
            swap_bytes((char*) &x, sizeof(Flt64));
            swap_bytes((char*) &y, sizeof(Flt64));
            swap_bytes((char*) &z, sizeof(Flt64));
        }

        bcopy((void*) &x, (void*) &buf[cur_offset], sizeof(Flt64));
        cur_offset += sizeof(Flt64);
        bcopy((void*) &y, (void*) &buf[cur_offset], sizeof(Flt64));
        cur_offset += sizeof(Flt64);
        bcopy((void*) &z, (void*) &buf[cur_offset], sizeof(Flt64));
        cur_offset += sizeof(Flt64);
        return cur_offset;
    }

    Int32 anyBuf::pushLittleEndianInt16(const Int16 iv)
    {
        Int16 v = iv;

        if (bigEndian())
        {
            swap_bytes((char*) &v, 2);
        }

        bcopy((void*) &v, (void*) &buf[cur_offset], 2);
        cur_offset += 2;
        return cur_offset;
    }


    Int32 anyBuf::pushLittleEndianUInt16(const UInt16 iv)
    {
        UInt16 v = iv;

        if (bigEndian())
        {
            swap_bytes((char*) &v, 2);
        }

        bcopy((void*) &v, (void*) &buf[cur_offset], 2);
        cur_offset += 2;
        return cur_offset;
    }
    */

    int32_t pushLittleEndianInt32(const int32_t iv)
    {
        int32_t v = iv;

        if (bigEndian())
        {
            swap_bytes((uint8_t*)&v, 4);
        }

        bcopy((void*)&v, (void*)&buf[cur_offset], 4);
        cur_offset += 4;
        return cur_offset;
    }
    /*

    Int32 anyBuf::pushLittleEndianS64(const s64 iv)
    {
        s64 v = iv;

        if (bigEndian())
        {
            swap_bytes((char*) &v, 8);
        }

        bcopy((void*) &v, (void*) &buf[cur_offset], 8);
        cur_offset += 8;
        return cur_offset;
    }

    */
    int32_t pushLittleEndianU64(const uint64_t iv)
    {
        uint64_t v = iv;

        if (bigEndian())
        {
            swap_bytes((uint8_t*)&v, 8);
        }

        bcopy((void*)&v, (void*)&buf[cur_offset], 8);
        cur_offset += 8;
        return cur_offset;
    }

    /*
    Int32 anyBuf::pushLittleEndianFlt32(const Flt32 iv)
    {
        Flt32 v = iv;

        if (bigEndian())
        {
            swap_bytes((char*) &v, 4);
        }

        bcopy((void*) &v, (void*) &buf[cur_offset], 4);
        cur_offset += 4;
        return cur_offset;
    }

    */
    int32_t pushLittleEndianFlt64(const double iv)
    {
        double v = iv;

        if (bigEndian())
        {
            swap_bytes((uint8_t*)&v, 8);
        }

        bcopy((void*)&v, (void*)&buf[cur_offset], 8);
        cur_offset += 8;
        return cur_offset;
    }

private:
    int32_t pop(int32_t bytes, void* result, bool swap)
    {
        bcopy((void*)&buf[cur_offset], result, bytes);
        cur_offset += bytes;

        if (swap)
        {
            swap_bytes((uint8_t*)result, bytes);
        }

        return 0;
    }

    void swap_bytes(uint8_t* data, int type)
    {
        char tmp;
        int i, j;

        for (i = 0, j = type - 1; i < j; i++, j--)
        {
            tmp = data[i];
            data[i] = data[j];
            data[j] = tmp;
        }
    }

    int32_t max_size;
    int32_t cur_size;
    int32_t cur_offset;
    uint8_t* buf;
};
#endif
