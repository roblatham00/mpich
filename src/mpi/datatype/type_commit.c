/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_Type_commit */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Type_commit = PMPI_Type_commit
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Type_commit  MPI_Type_commit
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Type_commit as PMPI_Type_commit
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Type_commit(MPI_Datatype * datatype) __attribute__ ((weak, alias("PMPI_Type_commit")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Type_commit
#define MPI_Type_commit PMPI_Type_commit

/*@
  MPIR_Type_commit

Input Parameters:
. datatype_p - pointer to MPI datatype

Output Parameters:

  Return Value:
  0 on success, -1 on failure.
@*/
int MPIR_Type_commit(MPI_Datatype * datatype_p)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_Datatype *datatype_ptr;

    MPIR_Assert(HANDLE_GET_KIND(*datatype_p) != HANDLE_KIND_BUILTIN);

    MPIR_Datatype_get_ptr(*datatype_p, datatype_ptr);

    if (datatype_ptr->is_committed == 0) {
        datatype_ptr->is_committed = 1;

#ifdef WITH_DAME
        /* Although we pass dataloop_size by pointer, it will not get set in
         * the function. It was relatively easy to calculate the size when
         * there were no significant optimizations of the datatype taking
         * place when constructing the dataloop. But now that we aggressively
         * optimize the representation, calculating the size becomes difficult.
         * So we run a second pass to calculate the actual size after all the
         * optimizations have been performed. I'll continue to pass the
         * parameter in case we later decide that we would rather not do
         * another pass over the datatype and calculate the size "on the fly" */
        MPIR_Dame_create(*datatype_p,
                         &datatype_ptr->dataloop,
                         &datatype_ptr->dataloop_size, &datatype_ptr->dataloop_depth);

        /* We need the dataloop size to compact it for the RMA code. */
        MPIR_Dame_calculate_size(datatype_ptr->dataloop, &datatype_ptr->dataloop_size, 0);

        /* After the call to Dataloop_create, the structs will be represented
         * without the additional space above the subtypes which DL_pack and
         * DL_unpack expect. While it's in this form, it's easier (and more
         * space-efficient) to compact it. The call to Dataloop_update() will
         * "fix" the dataloop to the form that DL_pack and DL_unpack expect.
         * (for details see the comments in Dataloop_update()) */
        MPIR_Dame_serialize(datatype_ptr->dataloop,
                            datatype_ptr->dataloop_size, &datatype_ptr->compact_dataloop);

        /* This will adjust the inner types. See comments above and in the
         * definition of Dataloop_update. This representation will be used for
         * local operations and pt2pt communication. The compact representation
         * will be used in RMA operations and converted to this representation
         * on the target */
        MPIR_Dame_update(datatype_ptr->dataloop, 0);
#if 0
        MPIR_Dame_print(datatype_ptr->dataloop);
        MPIR_Dame_print_compact(datatype_ptr->compact_dataloop);
#endif
#else
#ifdef MPID_NEEDS_DLOOP_ALL_BYTES
        /* If MPID implementation needs use to reduce everything to
         * a byte stream, do that. */
        MPIR_Dataloop_create(*datatype_p,
                             &datatype_ptr->dataloop,
                             &datatype_ptr->dataloop_size,
                             &datatype_ptr->dataloop_depth, MPIDU_DATALOOP_ALL_BYTES);
#else
        MPIR_Dataloop_create(*datatype_p,
                             &datatype_ptr->dataloop,
                             &datatype_ptr->dataloop_size,
                             &datatype_ptr->dataloop_depth, MPIR_DATALOOP_HOMOGENEOUS);
#endif

        /* create heterogeneous dataloop */
        MPIR_Dataloop_create(*datatype_p,
                             &datatype_ptr->hetero_dloop,
                             &datatype_ptr->hetero_dloop_size,
                             &datatype_ptr->hetero_dloop_depth, MPIR_DATALOOP_HETEROGENEOUS);

        MPL_DBG_MSG_D(MPIR_DBG_DATATYPE, TERSE, "# contig blocks = %d\n",
                      (int) datatype_ptr->max_contig_blocks);

#if 0
        MPII_Dataloop_dot_printf(datatype_ptr->dataloop, 0, 1);
#endif
#endif
#ifdef MPID_Type_commit_hook
        MPID_Type_commit_hook(datatype_ptr);
#endif /* MPID_Type_commit_hook */

    }
    return mpi_errno;
}

#undef FUNCNAME
#define FUNCNAME MPIR_Type_commit_impl
#undef FCNAME
#define FCNAME MPL_QUOTE(FUNCNAME)
int MPIR_Type_commit_impl(MPI_Datatype * datatype)
{
    int mpi_errno = MPI_SUCCESS;

    if (HANDLE_GET_KIND(*datatype) == HANDLE_KIND_BUILTIN)
        goto fn_exit;

    /* pair types stored as real types are a special case */
    if (*datatype == MPI_FLOAT_INT ||
        *datatype == MPI_DOUBLE_INT ||
        *datatype == MPI_LONG_INT || *datatype == MPI_SHORT_INT || *datatype == MPI_LONG_DOUBLE_INT)
        goto fn_exit;

    mpi_errno = MPIR_Type_commit(datatype);
    if (mpi_errno)
        MPIR_ERR_POP(mpi_errno);

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

#endif

#undef FUNCNAME
#define FUNCNAME MPI_Type_commit
#undef FCNAME
#define FCNAME MPL_QUOTE(FUNCNAME)
/*@
    MPI_Type_commit - Commits the datatype

Input Parameters:
. datatype - datatype (handle)

.N ThreadSafe

.N Fortran

.N Errors
.N MPI_SUCCESS
.N MPI_ERR_TYPE
@*/
int MPI_Type_commit(MPI_Datatype * datatype)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_FUNC_TERSE_STATE_DECL(MPID_STATE_MPI_TYPE_COMMIT);

    MPIR_ERRTEST_INITIALIZED_ORDIE();

    MPID_THREAD_CS_ENTER(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);
    MPIR_FUNC_TERSE_ENTER(MPID_STATE_MPI_TYPE_COMMIT);

    /* Validate parameters, especially handles needing to be converted */
#ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            MPIR_ERRTEST_ARGNULL(datatype, "datatype", mpi_errno);
            MPIR_ERRTEST_DATATYPE(*datatype, "datatype", mpi_errno);
        }
        MPID_END_ERROR_CHECKS;
    }
#endif

    /* Validate parameters and objects (post conversion) */
#ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            MPIR_Datatype *datatype_ptr = NULL;

            /* Convert MPI object handles to object pointers */
            MPIR_Datatype_get_ptr(*datatype, datatype_ptr);

            /* Validate datatype_ptr */
            MPIR_Datatype_valid_ptr(datatype_ptr, mpi_errno);
            if (mpi_errno)
                goto fn_fail;
        }
        MPID_END_ERROR_CHECKS;
    }
#endif /* HAVE_ERROR_CHECKING */

    /* ... body of routine ... */

    mpi_errno = MPIR_Type_commit_impl(datatype);
    if (mpi_errno)
        MPIR_ERR_POP(mpi_errno);

    /* ... end of body of routine ... */

  fn_exit:
    MPIR_FUNC_TERSE_EXIT(MPID_STATE_MPI_TYPE_COMMIT);
    MPID_THREAD_CS_EXIT(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);
    return mpi_errno;

  fn_fail:
    /* --BEGIN ERROR HANDLING-- */
#ifdef HAVE_ERROR_CHECKING
    {
        mpi_errno =
            MPIR_Err_create_code(mpi_errno, MPIR_ERR_RECOVERABLE, FCNAME, __LINE__, MPI_ERR_OTHER,
                                 "**mpi_type_commit", "**mpi_type_commit %p", datatype);
    }
#endif
    mpi_errno = MPIR_Err_return_comm(NULL, FCNAME, mpi_errno);
    goto fn_exit;
    /* --END ERROR HANDLING-- */
}
