#include <assert.h>

#include "adio.h"
#include "adioi.h"
#include "logfs_rtree.h"
#include "memstack.h"
#include "growvector.h"



/* struct uniquely describing rtree item */
typedef struct {
    rtree_range range;
    ADIO_Offset diskstart;
} logfs_rtree_item;


typedef struct {
    int loops;                  /* number of write stages we need */
    int globalloops;            /* global max number of loops needed */
    int coll;                   /* collective mode ? */
    void *cbdata;
    const logfs_rtree_flush_cb *cb;
    MPI_Comm comm;
    int bufsize;
    void *readbuf;
    void *writebuf;
    int done;                   /* flag to see if we're done */
    ADIO_Offset filesize;


    /* accept state */
    growvector_handle blocklens;        /* size of region */
    growvector_handle indices;  /* position in datafile */
    growvector_handle realpos;  /* position in real file */
    ADIO_Offset writesize;      /* amount of data to write so far */

} logfs_rtree_flush_state;

static inline void myqsort_swap(MPI_Aint numbers1[], MPI_Count numbers2[], int s1, int s2)
{
    MPI_Count tmp;
    MPI_Aint tmpa;
    tmpa = numbers1[s1];
    numbers1[s1] = numbers1[s2];
    numbers1[s2] = tmpa;
    tmp = numbers2[s1];
    numbers2[s1] = numbers2[s2];
    numbers2[s2] = tmp;
}

static inline void myqsort(MPI_Aint numbers1[], MPI_Count numbers2[], int left, int right)
{
    int l_hold = left;
    int r_hold = right;
    MPI_Aint pivot1 = numbers1[left];
    MPI_Count pivot2 = numbers2[left];
    while (left < right)
    {
        while ((numbers1[right] >= pivot1) && (left < right))
            right--;
        if (left != right)
            myqsort_swap(numbers1, numbers2, left, right);

        while ((numbers1[left] <= pivot1) && (left < right))
            left++;
        if (left != right) {
            myqsort_swap(numbers1, numbers2, left, right);
            right--;
        }
    }
    numbers1[left] = pivot1;
    numbers2[left] = pivot2;
    pivot1 = left;
    left = l_hold;
    right = r_hold;
    if (left < pivot1)
        myqsort(numbers1, numbers2, left, pivot1 - 1);
    if (right > pivot1)
        myqsort(numbers1, numbers2, pivot1 + 1, right);
}

/* takes the growvectors and build list sorted on indices so that the
 * resulting list can be used to construct the filetype for the datalog read
 */
static void logfs_rtree_readtypes(logfs_rtree_flush_state *state, MPI_Datatype *memtype,
                                  MPI_Datatype *filetype, int count, MPI_Aint *sortindices,
                                  MPI_Count *sortblocklens)
{
    ADIO_Offset *blocklens = growvector_get_null(state->blocklens);
    int i;
    ADIO_Offset *tmp1 = (ADIO_Offset *) growvector_get_null(state->indices);
    ADIO_Offset *tmp2 = (ADIO_Offset *) growvector_get_null(state->blocklens);

    /* TODO here: see if we can join readblocks */

    /* convert to ints */
    for (i = 0; i < count; ++i)
    {
        sortindices[i] = (MPI_Aint)(*tmp1++);
        sortblocklens[i] = (MPI_Count)(*tmp2++);
    }

    /* Moet ook voor originele volgorde de gesorteerde index weten om memtype
     * aan te maken */
   /* google translate: Must also know the sorted index for the original order
    * to create memtype */
    myqsort(sortindices, sortblocklens, 0, count - 1);
    ADIOI_Type_create_hindexed_x(count, sortblocklens, sortindices, MPI_BYTE, filetype);
    MPI_Type_commit(filetype);

    /* now reuse sortindices to create memtype */
    /* WAARSCHIJNLIJK VERKEERD!!! Dit doet geen herordering in geheugen! */
   /* google translate: PROBABLY WRONG This does not Reordered the memory! */
    sortindices[0] = 0;
    sortblocklens[0] = blocklens[0];
    for (i = 1; i < count; ++i) {
        sortindices[i] = sortindices[i - 1] + blocklens[i - 1];
        sortblocklens[i] = blocklens[i];
    }

    ADIOI_Type_create_hindexed_x(count, sortblocklens, sortindices, MPI_BYTE, memtype);
    MPI_Type_commit(memtype);


}

/* take the lists in blocklens, indices & realpos and try to initiate a
 * read operation; if the read is ready, try to start a write operation */
static void logfs_rtree_replay_startwrite(logfs_rtree_flush_state * state)
{
    MPI_Datatype readmemtype;
    MPI_Datatype readfiletype;
    MPI_Datatype writefiletype;
    int segmentcount = growvector_size(state->indices);

    /* temp buffers for ADIO_Offset -> int conversion */
    MPI_Aint *sortindices = (MPI_Aint *)ADIOI_Malloc(sizeof(MPI_Aint) * segmentcount);
    MPI_Count *sortblocklens = (MPI_Count *)ADIOI_Malloc(sizeof(MPI_Count) * segmentcount);

    void *buffer = state->readbuf;
    ADIO_Offset *tmp1;
    ADIO_Offset *tmp2;
    int i;

    /* create the read types:
     *   - a (sorted) list for reading from the datalog
     *   - a reordering memory type (putting data in order of the dataorder in
     *   the real file)
     *   - reuses buffers created above
     **/

    /* If we don't have anything to do, we don't have to read our logfile
     * (since that one is private to our CPU)
     * but we DO have to write if we are collective
     */
    if (segmentcount) {

        logfs_rtree_readtypes(state, &readmemtype, &readfiletype, segmentcount,
                              sortindices, sortblocklens);

        /* read */
        state->cb->readstart(buffer, readmemtype, readfiletype, state->cbdata);
        state->cb->readwait(state->cbdata);

        MPI_Type_free(&readmemtype);
        MPI_Type_free(&readfiletype);

        /*
         * if needed, we need to do 0-byte collective writes
         * (but not in this function) until everybody is finished
         */


        /* filetype for real file; need to convert to ints  */
        tmp1 = (ADIO_Offset *)growvector_get_null(state->blocklens);
        tmp2 = (ADIO_Offset *)growvector_get_null(state->realpos);
        for (i = 0; i < segmentcount; ++i)
        {
            sortblocklens[i] = (MPI_Count)(*tmp1++);
            sortindices[i] = (MPI_Aint)(*tmp2++);
        }

        ADIOI_Type_create_hindexed_x(segmentcount, sortblocklens, sortindices, MPI_BYTE, &writefiletype);
        MPI_Type_commit(&writefiletype);

    }
    else {
        /* create dummy write type */
        MPI_Type_contiguous(0, MPI_BYTE, &writefiletype);
        MPI_Type_commit(&writefiletype);
    }

    /* write */
    state->cb->writestart(buffer, writefiletype, state->writesize, state->cbdata);
    state->cb->writewait(state->cbdata);
    state->writesize = 0;
    MPI_Type_free(&writefiletype);

    growvector_clear(state->blocklens);
    growvector_clear(state->indices);
    growvector_clear(state->realpos);

    /* free conversion buffers */
    ADIOI_Free(sortindices);
    ADIOI_Free(sortblocklens);
    
    /* if we are collective, do an alreduce in the end to find the maximum
     * filesize and to see if everybody is done... */
}

/* try to add region; Return remaining bytes if full */
/* call with rangestart == rangestop == 0 for indicating the end of data */
/* probably need to add extra vars here: alldone and filesize */
static ADIO_Offset logfs_rtree_flush_add(logfs_rtree_flush_state * state,
                                         ADIO_Offset rangestart, ADIO_Offset rangestop,
                                         ADIO_Offset fileofs)
{
    ADIO_Offset leftover = 0;
    ADIO_Offset thiswrite = rangestop - rangestart;
    int forcewrite = 0;

    /* accept data until we hit the buffer size */
    /* (maybe split access) */
    /* reduce size if needed */
    if (thiswrite + state->writesize > state->bufsize) {
        thiswrite = (state->bufsize - state->writesize);
        leftover = rangestop - rangestart - thiswrite;
        rangestop = rangestart + thiswrite;

        /* could do something smart here and decide it is not worth it to split
         * up a small region into even smaller ones (e.g. splitting a 8 byte
         * write into 2 regions...) */
    }

    if (thiswrite) {
        /* add write to lists */
        growvector_pushback(state->blocklens, &thiswrite, sizeof(thiswrite));
        growvector_pushback(state->indices, &fileofs, sizeof(fileofs));
        growvector_pushback(state->realpos, &rangestart, sizeof(fileofs));
        state->writesize += thiswrite;
    }
    else {
        /* thiswrite == 0, meaning the user passed rangestart == rangestop
         * -> no more data, force write */
        forcewrite = 1;
    }

    /* progress reads & writes */

    if (state->writesize < state->bufsize && !forcewrite)
        return leftover;

    /* OK, databuffer is full or there will be no more data
     * Do write */
    logfs_rtree_replay_startwrite(state);
   /* filled up this buffer: will issue a collective write operation with
    * real data, so we have one less no-op write to issue */
   state->loops--;

    return leftover;
}


/* accept a region from the rtree; Add the data to the list of data
 * to be written; If the buffer size is full, redo with the rest of the region
 * */
static int logfs_rtree_flush_accept(const rtree_range * range, ADIO_Offset * fileofs, void *extra)
{
    logfs_rtree_flush_state *state = (logfs_rtree_flush_state *) extra;

   /* all the bytes described by this node of the tree */
    ADIO_Offset todo = range->stop - range->start;
   /* bytes processed by slave file system */
    ADIO_Offset done = 0;
   /* bytes processed in a specific flush_add call. Large ranges will require
    * multiple rounds */
   ADIO_Offset nbytes;

   /* this loop is semantically similar to that in ad_write.c, where we need to
    * write N bytes but can only issue a smaller number at a time */
   while (todo)
   {
       /* a small range that does not fill up the buffer could be handled in
	* one shot.  In fact, it might take multiple tree nodes to fill up the
	* buffer.  large ranges, though, will require multiple rounds */
      todo = logfs_rtree_flush_add (state, range->start + done,
	    range->stop, (*fileofs) + done);
      nbytes = range->stop - (range->start + done) - todo;
      done += nbytes;
    }
   /* after loop exits, there is likely (unless perfectly multiple of buffer
    * size) a partial write.  calling code will send an "end of data" call,
    * which will flush everything outstanding . */

    return 1;
}


void logfs_rtree_flush(logfs_rtree * tree, int bufsize,
                       const logfs_rtree_flush_cb * cb, void *cbdata, int coll,
                       ADIO_Offset * filesize, MPI_Comm comm)
{
    logfs_rtree_flush_state state;

   if (!tree->rtree) return;

    state.coll = coll;
    state.cbdata = cbdata;
    state.cb = cb;
    state.comm = comm;
    state.bufsize = bufsize;
    state.readbuf = ADIOI_Malloc(bufsize);
    state.writebuf = ADIOI_Malloc(bufsize);
    state.globalloops = 0;
    state.done = 0;
    state.filesize = *filesize;
    state.writesize = 0;

    /* allocate arrays; preallocate 1k for each */
    state.blocklens = growvector_create(sizeof(ADIO_Offset), (1024 / 8));
    state.indices = growvector_create(sizeof(ADIO_Offset), (1024 / 8));
    state.realpos = growvector_create(sizeof(ADIO_Offset), (1024 / 8));

    if (state.coll) {
        int buf = bufsize;
#ifndef NDEBUG
        MPI_Bcast(&buf, 1, MPI_INT, 0, comm);
        assert(buf == bufsize);
#endif
    }

    assert(bufsize);

    /* calculate  number of write stages */
    state.loops = (tree->rangesize + bufsize - 1) / bufsize;

    if (coll) {
        MPI_Allreduce(&state.loops, &state.globalloops, 1, MPI_INT, MPI_MAX, comm);
        if (state.globalloops > state.loops)
            state.loops = state.globalloops;
    }

    /* Notify caller that we are going to start  */
    state.cb->start(state.cbdata, state.coll);

    /* for now go over tree in order */
    rtree_walk(tree->rtree, logfs_rtree_flush_accept, &state);

    /* indicate end of data by writing empty region  */
    logfs_rtree_flush_add(&state, 0, 0, 0);

   if (state.coll) {
       while (state.loops > 0) {
	   /* Possible for I/o workloaod to be imbalanced across processes.
	    * Issue zero-byte writes until other processes have finished their
	    * collective writes */
	   state.cb->writestart(NULL, MPI_BYTE, 0, state.cbdata);
	   state.cb->writewait(state.cbdata);
	   state.loops--;
       }
   }


    /* shut down callbacks */
    state.cb->stop(state.cbdata);

    /* cleanup */
    growvector_free(&state.blocklens);
    growvector_free(&state.indices);
    growvector_free(&state.realpos);

    /* === */
    ADIOI_Free(state.readbuf);
    ADIOI_Free(state.writebuf);

    /* update filesize */
    if (coll)
        MPI_Allreduce(&state.filesize, filesize, 1, ADIO_OFFSET, MPI_MAX, comm);
    else
        *filesize = state.filesize;
}



/* rtree callback that adds the item to the list */
static int rtree_add_item(const rtree_range * range, ADIO_Offset * diskstart, void *extra)
{
    memstack_handle list = (memstack_handle) extra;
    logfs_rtree_item *newitem = (logfs_rtree_item *)
        memstack_push(list);

    /* copy item into memstack */
    newitem->diskstart = *diskstart;
    newitem->range = *range;

    return 1;
}

static inline void logfs_rtree_updatesize(ADIO_Offset * datasize, ADIO_Offset delta)
{
    if (!datasize)
        return;

    if (delta < 0) {
        assert(*datasize >= -delta);
    }
    *datasize += delta;
}

/* TODO: merge adjacent regions if possible!
 *  Could be done by expanding the region by 1 on the left side (cannot
 *  normally merge on the right side since the datalog offsets only increase)
 *  If we find something, check if the datalog regions also match; if so,
 *  expand existing region */
void logfs_rtree_addsplit(logfs_rtree * tree, ADIO_Offset start, ADIO_Offset
                          stop, ADIO_Offset diskstart)
{
    rtree_range newrange;
    memstack_handle list;
    int listsize;
    rtree_handle rtree = tree->rtree;
    ADIO_Offset *datasize = &tree->rangesize;

    newrange.start = start;
    newrange.stop = stop;

    list = memstack_create(sizeof(logfs_rtree_item));

#if 0
    if (newrange.start > 0)
        --newrange.start;
#endif

    /* find regions in tree overlapping our new region */
    rtree_overlap(rtree, &newrange, rtree_add_item, list);

    /* process the search results from the list, removing the region from the
     * tree if it is completely inside the new region, and splitting it
     * otherwise */
    listsize = memstack_getsize(list);
    while (listsize--) {
        logfs_rtree_item *cur = (logfs_rtree_item *) memstack_pop(list);
        assert(cur);

        /* In any case the region needs to be removed */
        rtree_remove(rtree, &cur->range, 0);
        logfs_rtree_updatesize(datasize, -(cur->range.stop - cur->range.start));

        /* find out if it is completely inside; if so skip it  */
        if (cur->range.start >= newrange.start && cur->range.stop <= newrange.stop)
            continue;

        /* need to split the region */
        /* calculate split on the left */
        if (cur->range.start < newrange.start) {
            rtree_range new;
            new.start = cur->range.start;
            new.stop = newrange.start;
            assert(new.start < new.stop);

            /* add left part */
            rtree_add(rtree, &new, cur->diskstart);
            logfs_rtree_updatesize(datasize, new.stop - new.start);
        }
        if (cur->range.stop > newrange.stop) {
            rtree_range new;
            new.start = newrange.stop;
            new.stop = cur->range.stop;
            assert(new.start < new.stop);

            /* add right part, shifting diskstart
             * (except when the diskstart is invalid)*/
            rtree_add(rtree, &new,
                      (cur->diskstart = ADIO_OFFSET_INVALID ?
                       ADIO_OFFSET_INVALID : cur->diskstart + (new.start - cur->range.start)));
            logfs_rtree_updatesize(datasize, new.stop - new.start);
        }

        /*
         *  IDEA: could also record 'free' items in logfile; could be reused but
         *  will cause disk seeking!
         */
    }

    /* add the region */
    rtree_add(rtree, &newrange, diskstart);
    logfs_rtree_updatesize(datasize, newrange.stop - newrange.start);
   memstack_free(&list);
}
