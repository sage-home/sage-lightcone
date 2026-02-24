#ifndef LIBHPC_MPI_STUB_HH
#define LIBHPC_MPI_STUB_HH

// Check if a real MPI header has been included (e.g. by HDF5)
// OR if using parallel HDF5 in serial build mode
#if defined(MPI_VERSION) || defined(MPI_H) || defined(OMPI_MPI_H) || \
    defined(USE_PARALLEL_HDF5_IN_SERIAL)

#if defined(USE_PARALLEL_HDF5_IN_SERIAL) && !defined(MPI_VERSION)
#include <mpi.h>
#endif

// If real MPI types are present, we generally rely on them.
// If we are not linking to the real MPI library, we might have linker errors.
// However, typically if one includes mpi.h, one links -lmpi.
// For the purpose of "removing mpi", ensure we use the serial HDF5 which won't
// include mpi.h.

#else

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <vector>

// Types
typedef int MPI_Comm;
typedef int MPI_Info;
typedef int MPI_Datatype;
typedef int MPI_Op;
typedef int MPI_Request;
typedef int MPI_Group;

struct MPI_Status
{
    int MPI_SOURCE;
    int MPI_TAG;
    int MPI_ERROR;
};

// Constants
#define MPI_COMM_WORLD 0
#define MPI_COMM_SELF 1
#define MPI_COMM_NULL 2
#define MPI_INFO_NULL 0

#define MPI_SUCCESS 0
#define MPI_ERR_OTHER 1
#define MPI_MAX_PROCESSOR_NAME 256
#ifndef MPI_STATUSES_IGNORE
#define MPI_STATUSES_IGNORE ((MPI_Status*)0)
#endif
#ifndef MPI_STATUS_IGNORE
#define MPI_STATUS_IGNORE ((MPI_Status*)0)
#endif

#define MPI_INT 1
#define MPI_LONG 2
#define MPI_LONG_LONG 3
#define MPI_UNSIGNED 4
#define MPI_UNSIGNED_LONG 5
#define MPI_UNSIGNED_LONG_LONG 6
#define MPI_FLOAT 7
#define MPI_DOUBLE 8
#define MPI_CHAR 9
#define MPI_BYTE 10
#define MPI_C_BOOL 11

#define MPI_SUM 1
#define MPI_MAX 2
#define MPI_MIN 3
#define MPI_PROD 4

#define MPI_ANY_SOURCE -1
#define MPI_ANY_TAG -1

#define MPI_DATATYPE_NULL 0

// Use ifndef/def to prevent redefinition if these macros leaked from somewhere
// else
#ifndef MPI_MODE_CREATE
#define MPI_MODE_CREATE 1
#endif
#ifndef MPI_MODE_WRONLY
#define MPI_MODE_WRONLY 2
#endif
#ifndef MPI_MODE_RDONLY
#define MPI_MODE_RDONLY 4
#endif

#define MPI_REQUEST_NULL 0
#define MPI_IN_PLACE ((void*)-1)
#define MPI_UNDEFINED -1

// Type Helpers (forward)
inline int MPI_Type_size(MPI_Datatype datatype, int* size);

// Function Stubs
inline int MPI_Init(int* argc, char*** argv) { return MPI_SUCCESS; }
inline int MPI_Finalize() { return MPI_SUCCESS; }
inline int MPI_Initialized(int* flag)
{
    *flag = 1;
    return MPI_SUCCESS;
}

inline int MPI_Comm_rank(MPI_Comm comm, int* rank)
{
    *rank = 0;
    return MPI_SUCCESS;
}

inline int MPI_Comm_size(MPI_Comm comm, int* size)
{
    *size = 1;
    return MPI_SUCCESS;
}

inline int MPI_Barrier(MPI_Comm comm) { return MPI_SUCCESS; }

inline int MPI_Comm_dup(MPI_Comm comm, MPI_Comm* newcomm)
{
    *newcomm = comm;
    return MPI_SUCCESS;
}

inline int MPI_Comm_free(MPI_Comm* comm)
{
    *comm = MPI_COMM_NULL;
    return MPI_SUCCESS;
}

inline int MPI_Comm_group(MPI_Comm comm, MPI_Group* group)
{
    *group = 0;
    return MPI_SUCCESS;
}

inline int MPI_Group_free(MPI_Group* group)
{
    *group = 0;
    return MPI_SUCCESS;
}

inline int MPI_Group_translate_ranks(MPI_Group group1, int n, const int* ranks1, MPI_Group group2,
                                     int* ranks2)
{
    for (int i = 0; i < n; ++i)
        ranks2[i] = ranks1[i];
    return MPI_SUCCESS;
}

inline int MPI_Group_excl(MPI_Group group, int n, const int* ranks, MPI_Group* newgroup)
{
    *newgroup = 0;
    return MPI_SUCCESS;
}

inline int MPI_Group_range_incl(MPI_Group group, int n, int (*ranges)[3], MPI_Group* newgroup)
{
    *newgroup = 0;
    return MPI_SUCCESS;
}

inline int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm* newcomm)
{
    *newcomm = comm;
    return MPI_SUCCESS;
}

inline int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm* newcomm)
{
    *newcomm = comm;
    return MPI_SUCCESS;
}

inline int MPI_Sendrecv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest,
                        int sendtag, void* recvbuf, int recvcount, MPI_Datatype recvtype,
                        int source, int recvtag, MPI_Comm comm, MPI_Status* status)
{
    if (dest == 0 && source == 0)
    {
        int size = 0;
        MPI_Type_size(sendtype, &size);
        std::memcpy(recvbuf, sendbuf, sendcount * size);
    }
    return MPI_SUCCESS;
}

inline int MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status* status)
{
    if (status)
    {
        status->MPI_SOURCE = source;
        status->MPI_TAG = tag;
        status->MPI_ERROR = MPI_SUCCESS;
    }
    return MPI_SUCCESS;
}

inline int MPI_Abort(MPI_Comm comm, int errorcode)
{
    std::cerr << "MPI_Abort called with error " << errorcode << std::endl;
    std::exit(errorcode);
    return MPI_SUCCESS;
}

// Type Helpers
inline int MPI_Type_size(MPI_Datatype datatype, int* size)
{
    if (datatype == MPI_INT)
        *size = sizeof(int);
    else if (datatype == MPI_LONG)
        *size = sizeof(long);
    else if (datatype == MPI_LONG_LONG)
        *size = sizeof(long long);
    else if (datatype == MPI_UNSIGNED)
        *size = sizeof(unsigned);
    else if (datatype == MPI_UNSIGNED_LONG)
        *size = sizeof(unsigned long);
    else if (datatype == MPI_UNSIGNED_LONG_LONG)
        *size = sizeof(unsigned long long);
    else if (datatype == MPI_FLOAT)
        *size = sizeof(float);
    else if (datatype == MPI_DOUBLE)
        *size = sizeof(double);
    else if (datatype == MPI_CHAR)
        *size = sizeof(char);
    else if (datatype == MPI_BYTE)
        *size = 1;
    else if (datatype == MPI_C_BOOL)
        *size = sizeof(bool);
    else
        *size = 1;
    return MPI_SUCCESS;
}

inline int MPI_Type_commit(MPI_Datatype* datatype) { return MPI_SUCCESS; }
inline int MPI_Type_free(MPI_Datatype* datatype) { return MPI_SUCCESS; }
inline int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype* newtype)
{
    *newtype = oldtype; // Simplification
    return MPI_SUCCESS;
}
inline int MPI_Type_indexed(int count, const int* array_of_blocklengths,
                            const int* array_of_displacements, MPI_Datatype oldtype,
                            MPI_Datatype* newtype)
{
    *newtype = oldtype;
    return MPI_SUCCESS;
}

inline int MPI_Type_create_indexed_block(int count, int blocklength,
                                         const int* array_of_displacements, MPI_Datatype oldtype,
                                         MPI_Datatype* newtype)
{
    *newtype = oldtype;
    return MPI_SUCCESS;
}

// Collective Operations
// For serial execution, many collectives are no-ops or simple copies
inline int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
{
    return MPI_SUCCESS;
}

inline int MPI_Allreduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
                         MPI_Op op, MPI_Comm comm)
{
    if (sendbuf == MPI_IN_PLACE)
        return MPI_SUCCESS;

    int type_size;
    MPI_Type_size(datatype, &type_size);
    // In a single process, the result is just the input
    if (sendbuf != recvbuf)
    {
        if (recvbuf && sendbuf)
        {
            std::memcpy(recvbuf, sendbuf, count * type_size);
        }
    }
    return MPI_SUCCESS;
}

inline int MPI_Reduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
                      MPI_Op op, int root, MPI_Comm comm)
{
    return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}

inline int MPI_Scan(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                    MPI_Comm comm)
{
    return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}

inline int MPI_Exscan(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
                      MPI_Op op, MPI_Comm comm)
{
    // Exscan on 1 process is undefined/unused
    return MPI_SUCCESS;
}

inline int MPI_Gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
                      int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    if (sendbuf == MPI_IN_PLACE)
        return MPI_SUCCESS;

    // Assume sendcount/recvcount are compatible items
    int type_size;
    MPI_Type_size(sendtype, &type_size); // Assuming sendtype == recvtype

    if (sendbuf != recvbuf)
    {
        if (recvbuf && sendbuf)
        {
            std::memcpy(recvbuf, sendbuf, sendcount * type_size);
        }
    }
    return MPI_SUCCESS;
}

inline int MPI_Allgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
                         int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
    return MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, 0, comm);
}

inline int MPI_Scatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
                       int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    // Inverse of Gather, but for 1 proc it's the same copy
    return MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
}

inline int MPI_Allgatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
                          const int* recvcounts, const int* displs, MPI_Datatype recvtype,
                          MPI_Comm comm)
{
    // Simplified stub
    if (recvcounts && recvcounts[0] > 0)
        return MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcounts[0], recvtype, comm);
    return MPI_SUCCESS;
}

inline int MPI_Gatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
                       const int* recvcounts, const int* displs, MPI_Datatype recvtype, int root,
                       MPI_Comm comm)
{
    if (recvcounts && recvcounts[0] > 0)
        return MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcounts[0], recvtype, comm);
    return MPI_SUCCESS;
}

// Point to Point
inline int MPI_Send(const void* buf, int count, MPI_Datatype datatype, int dest, int tag,
                    MPI_Comm comm)
{
    return MPI_SUCCESS;
}
inline int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                    MPI_Status* status)
{
    return MPI_SUCCESS;
}
inline int MPI_Wait(MPI_Request* request, MPI_Status* status) { return MPI_SUCCESS; }
inline int MPI_Test(MPI_Request* request, int* flag, MPI_Status* status)
{
    *flag = 1;
    return MPI_SUCCESS;
}
inline int MPI_Cancel(MPI_Request* request) { return MPI_SUCCESS; }
inline int MPI_Request_free(MPI_Request* request) { return MPI_SUCCESS; }

// Async point-to-point
inline int MPI_Isend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag,
                     MPI_Comm comm, MPI_Request* request)
{
    *request = 0;
    return MPI_SUCCESS;
}
inline int MPI_Issend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag,
                      MPI_Comm comm, MPI_Request* request)
{
    *request = 0;
    return MPI_SUCCESS;
}
inline int MPI_Irecv(void* buf, int count, MPI_Datatype datatype, int source, int tag,
                     MPI_Comm comm, MPI_Request* request)
{
    *request = 0;
    return MPI_SUCCESS;
}

// Multiple request operations
inline int MPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[])
{
    return MPI_SUCCESS;
}
inline int MPI_Testall(int count, MPI_Request array_of_requests[], int* flag,
                       MPI_Status array_of_statuses[])
{
    *flag = 1;
    return MPI_SUCCESS;
}

// Processor name
inline int MPI_Get_processor_name(char* name, int* resultlen)
{
    const char* hostname = "localhost";
    strncpy(name, hostname, MPI_MAX_PROCESSOR_NAME - 1);
    name[MPI_MAX_PROCESSOR_NAME - 1] = '\0';
    *resultlen = strlen(name);
    return MPI_SUCCESS;
}

// Info
inline int MPI_Info_create(MPI_Info* info)
{
    *info = MPI_INFO_NULL;
    return MPI_SUCCESS;
}
inline int MPI_Info_set(MPI_Info info, const char* key, const char* value) { return MPI_SUCCESS; }
inline int MPI_Info_free(MPI_Info* info) { return MPI_SUCCESS; }

// File I/O Definitions (dummy types if not defined)
typedef int MPI_File;
#ifndef MPI_FILE_NULL
#define MPI_FILE_NULL 0
#endif

// Some HDF5 codes might use these
inline int MPI_File_open(MPI_Comm comm, const char* filename, int amode, MPI_Info info,
                         MPI_File* fh)
{
    *fh = 0;
    return MPI_SUCCESS;
}
inline int MPI_File_close(MPI_File* fh) { return MPI_SUCCESS; }
inline int MPI_File_set_size(MPI_File fh, int size) { return MPI_SUCCESS; }

// Non-blocking probe
inline int MPI_Iprobe(int source, int tag, MPI_Comm comm, int* flag, MPI_Status* status)
{
    *flag = 0; // No message available in serial mode
    return MPI_SUCCESS;
}

// All-to-all with varying types (no-op in serial)
inline int MPI_Alltoallw(const void* sendbuf, const int sendcounts[], const int sdispls[],
                         const MPI_Datatype sendtypes[], void* recvbuf, const int recvcounts[],
                         const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm)
{
    // In serial mode, just copy sendbuf to recvbuf if same size
    return MPI_SUCCESS;
}

// Timing function
inline double MPI_Wtime()
{
    // Return wall clock time in seconds
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

#endif // MPI_VERSION check
#endif // LIBHPC_MPI_STUB_HH
