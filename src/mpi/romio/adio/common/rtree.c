#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#ifdef USE_EFENCE
#include <efence.h>
#endif

#include <assert.h>
#include <string.h>

#include "rtree.h"


#define RTREE_CHILD_MIN 1        /* minimum children in node */
#define RTREE_CHILD_MAX 4        /* maximum number of children */ 

#ifndef FALSE
#define FALSE 0 
#endif

#ifndef TRUE
#define TRUE 1
#endif

/*
 * TODO:
 *    - Make MIN/MAX child configurable at tree compilation time?
 */
struct rtree_node; 

struct rtree_entry
{
   rtree_range 		range; 
   RTREE_DATA_TYPE 	data;  
} ; 

typedef struct rtree_entry rtree_entry; 

struct rtree_node
{
   rtree_range          range;
   struct rtree_node  * parent; 
   union 
   {
   	struct rtree_node *  child[RTREE_CHILD_MAX];
	struct rtree_entry * entry[RTREE_CHILD_MAX]; 
   }; 
} ; 

typedef struct rtree_node rtree_node; 


struct rtree
{
   struct rtree_node *  root;      /* pointer to root node */ 
   int                  depth;     /* index of lowest level */
   int                  count;      /* number of data items in the tree */

   rtree_callback       freefunc;  /* called when removing an item from the tree */
   void * 		freedata;  /* extra pointer passed to the free function */ 

   rtree_callback_split splitfunc; /* pointer to split function */
   void *               splitdata; /* extra data */ 
}; 

typedef struct rtree rtree; 

typedef int (*rtree_callback_visit) (rtree_const_handle tree, const rtree_node * node, 
      int depth, void * extra); 

 /*======== iterator functions ===========*/

struct rtree_iterator
{
   int          createdepth; 
   int *        childnum;    
   rtree_node * node;
   rtree_handle tree; 
   int          depth; 
}; 

typedef struct rtree_iterator rtree_iterator; 

/****************************************************************************
 * Overlap helper functions                                                 *
 ****************************************************************************/

static inline RTREE_RANGE_TYPE rtree_range_type_min 
		(RTREE_RANGE_TYPE a, RTREE_RANGE_TYPE b); 

static inline RTREE_RANGE_TYPE rtree_range_type_max 
		(RTREE_RANGE_TYPE a, RTREE_RANGE_TYPE b); 

/* Return true if the ranges intersect */ 
static inline int rtree_range_has_overlap (const rtree_range * range1,
      			const rtree_range * range2);

/* If there is overlap, return true and return common range*/  
static inline int rtree_range_shared (const rtree_range * range1, 
      const rtree_range * range2, rtree_range * dest); 

/* Calculate the smallest range containing both rectangles */ 
static inline void rtree_range_extent (const rtree_range * range1, 
      const rtree_range * range2, rtree_range * dest);

/* Extend range by adding range2; Store result in range1 */ 
static inline void rtree_range_extend (rtree_range * range1, 
      const rtree_range * range2); 

/* Order rectangles, range1 first */ 
static inline void rtree_range_order (const rtree_range * * range1, 
      const rtree_range * * range2); 

/* return size of the range */
static inline RTREE_RANGE_TYPE rtree_range_size (const rtree_range * r); 

/* calculate how much the range size increases when adding the r2 to r1 */ 
static inline RTREE_RANGE_TYPE rtree_range_calc_extension (const rtree_range * r1, 
      const rtree_range * r2); 

/* Check if the range is empty */ 
static inline int rtree_range_empty (const rtree_range * r); 

/* check if the ranges are equal */
static inline int rtree_range_equals (const rtree_range * r1, const rtree_range * r2);

/* check if the first range contains all of the second range */
static inline int rtree_range_contains (const rtree_range * r1, const rtree_range * r2);

/****************************************************************************
 * Memory functions                                                         *
 ****************************************************************************/

static inline void  * rtree_mem_alloc_raw (int bytes)
{
   void * ptr = malloc (bytes); 
#ifdef RTREE_INIT_MEM
   if (ptr)
   	memset (ptr, 0, bytes);  
#endif
   return ptr; 
}

static inline void rtree_mem_free_raw (void  * ptr, int bytes)
{
#ifdef RTREE_CLEAR_MEM
   memset (ptr, 0, bytes); 
#endif
   free (ptr); 
}

static inline rtree_iterator * rtree_mem_alloc_iterator ()
{
   rtree_iterator * i = rtree_mem_alloc_raw (sizeof (rtree_iterator)); 
   return i; 
}

static inline void rtree_mem_free_iterator (rtree_iterator ** iter)
{
   assert (iter); 
   rtree_mem_free_raw (*iter, sizeof(rtree_iterator)); 
   *iter=0; 
}

static inline rtree_node * rtree_mem_alloc_node ()
{
   int i; 
   rtree_node * n = rtree_mem_alloc_raw (sizeof (rtree_node));

   n->range.start = n->range.stop = 0; 
   n->parent = 0; 
   for (i=0; i<RTREE_CHILD_MAX; ++i)
      n->child[i] = 0; 
   return n; 
}

static inline void rtree_mem_free_node (rtree_node ** node)
{
   assert (node); 
   rtree_mem_free_raw (*node, sizeof(rtree_node)); 
   *node=0; 
}

static inline  void rtree_mem_free_tree (rtree ** tree)
{
   assert (tree); 
   rtree_mem_free_raw (*tree, sizeof(rtree)); 
   *tree = 0; 
}

static inline rtree * rtree_mem_alloc_tree ()
{
   return rtree_mem_alloc_raw (sizeof (rtree)); 
}

static inline rtree_entry * rtree_mem_alloc_entry 
     (const rtree_range * range, RTREE_DATA_TYPE data)
{
   rtree_entry * new = rtree_mem_alloc_raw (sizeof (rtree_entry)); 
   new->range = *range; 
   new->data = data; 
   return new; 
}

static inline void rtree_mem_free_entry (rtree_entry ** entry)
{
   assert (entry); 
   rtree_mem_free_raw (*entry, sizeof(rtree_entry)); 
   *entry = 0; 
}




/****************************************************************************
 * Public rtree functions                                                   *
 ****************************************************************************/
static int rtree_node_compare_entry (const void * n1, const void * n2)
{
   if (*(const rtree_entry **)n1 == 0)
      return ( *(const rtree_node **) n2 == 0 ? 0 : -1); 
   
   if (*(const rtree_entry **)n2 == 0)
      /* we know that n1 cannot be 0 here */
      return 1; 

   /*  n1 and n2 != 0 */
   RTREE_RANGE_TYPE s1 = (*(const rtree_entry **)n1)->range.start; 
   RTREE_RANGE_TYPE s2 = (*(const rtree_entry **)n2)->range.start; 
   if (s1 < s2)
      return -1; 
   if (s1 > s2)
      return 1; 
   return 0; 
}


static int rtree_node_compare_node (const void * n1, const void * n2)
{
   if (*(const rtree_node **)n1 == 0)
      return ( *(const rtree_node **) n2 == 0 ? 0 : -1); 
   
   if (*(const rtree_node **)n2 == 0)
      /* we know that n1 cannot be 0 here */
      return 1; 

   /*  n1 and n2 != 0 */
   RTREE_RANGE_TYPE s1 = (*(const rtree_node **)n1)->range.start; 
   RTREE_RANGE_TYPE s2 = (*(const rtree_node **)n2)->range.start; 
   if (s1 < s2)
      return -1; 
   if (s1 > s2)
      return 1; 
   return 0; 
}

/* sort the child/entry pointers in a node so that iteraters will traverse
 * the tree in order and so intra-node searches can be made faster */
static void rtree_node_sort (rtree_node * node, int count, int leaf)
{
   /* probably very bad for already sorted sequences! Consider other
    * sorting algorithms */
   qsort ((leaf ? (void *) &node->child[0] : (void *) &node->entry[0]), 
         count, 
         sizeof (rtree_node *),
         (leaf ? rtree_node_compare_node : rtree_node_compare_entry)); 

}

static void rtree_find_internal (const rtree_node * node, 
      const rtree_range * range, int depth, int treedepth, 
      const rtree_node ** nodeptr, int * entrynum)
{
   int i; 
 
   assert (node); 
   assert (nodeptr); 

   if (depth != treedepth)
   {
      for (i=0; i<RTREE_CHILD_MAX; ++i)
      {
         if (!node->child[i])
            break; 

         if (rtree_range_contains (&node->child[i]->range, range))
         {
	     rtree_find_internal (node->child[i], 
                  range, depth+1, treedepth, nodeptr, entrynum);

	    /* if this subtree contains the item, return
	     * otherwise keep looking since another child could also contain
	     * the item */
	    if (*nodeptr) 
	       return; 
	 }
      }
   }
      else
   {
      for (i=0; i<RTREE_CHILD_MAX; ++i)
      {
         if (!node->entry[i])
            break; 

         if (rtree_range_equals(range, &node->entry[i]->range))
         {
            /* found item */ 
            *nodeptr = node;
            *entrynum = i; 
            return;
         }
      }
   }

   /* didn't find the range */ 
   return; 
}

RTREE_DATA_TYPE * rtree_find (rtree_const_handle tree, 
      const rtree_range * range)
{
   const rtree_node * nodeptr = 0;
   int entrynum = 0; 

   assert (tree); 
   rtree_find_internal (tree->root, range, 0, tree->depth, 
         &nodeptr, &entrynum); 

   if (nodeptr)
      return &nodeptr->entry[entrynum]->data; 
        else
      return 0; 
}


static int rtree_overlap_internal (const rtree_node * node, const rtree_range * range, 
      rtree_callback callback, int depth, int treedepth, void * extra)
{
   int i;

   assert (node); 

   if (depth == treedepth)
   {
      /* At lowest level; We have pointers to entries */
      for (i=0; i<RTREE_CHILD_MAX; ++i)
      {
	 if (!node->entry[i])
	 {
	    assert ((depth == 0 && i>=0) ||
		    (depth > 0 && i>=1)); 
	    break;
	 }

	 if (rtree_range_has_overlap (range, &node->entry[i]->range))
	 {
	    /* Found qualifying entry */ 
	    if (!callback (&node->entry[i]->range, &node->entry[i]->data, 
		     		extra))
	       return FALSE; 
	 }
      }
      return TRUE; 
   }

   /* Not at lowest level yet; we have child nodes */ 
   for (i=0; i<RTREE_CHILD_MAX; ++i)
   {
      if (!node->child[i])
      {
	 assert (i>0); 
	 break; 
      }
      if (rtree_range_has_overlap (range, &node->child[i]->range))
      {
	 /* search subtree */ 
	 if (!rtree_overlap_internal (node->child[i], range, callback, 
		  depth+1, treedepth, extra))
	    return FALSE; 
      }
   }
   return TRUE; 
}

/*
 * Check all children, ignore subtrees that do not overlap
 * Return false if the callback returned false
 * otherwise return true
 */
int rtree_overlap (rtree_const_handle tree, const rtree_range * range, 
                rtree_callback callback, void * extra)
{
   assert (tree->root); 

   return rtree_overlap_internal (tree->root, range, callback, 0, 
	 	tree->depth, extra); 
}


void rtree_free_node (rtree_const_handle tree, rtree_node * node, int depth, int treedepth)
{
   int i; 
   if (depth == treedepth)
   {
      /* node contains pointers to entries */
      for (i=0; i<RTREE_CHILD_MAX; ++i)
      {
	 if (!node->entry[i])
	    break; 

	 if (tree->freefunc)
	 {
	    tree->freefunc (&node->entry[i]->range, 
		  &node->entry[i]->data, tree->freedata);
	 }
	 rtree_mem_free_entry (&node->entry[i]);
      }
   }
   else
   {
      /* free all child nodes, then free our own node */
      for (i=0; i<RTREE_CHILD_MAX; ++i)
      {
	 if (!node->child[i])
	    break;
	 rtree_free_node (tree, node->child[i], depth+1, treedepth); 
      }
   }
   rtree_mem_free_node (&node); 
}

void rtree_free (rtree_handle * rtree)
{
    if (!rtree) return;
    if (*rtree == NULL) return;
    if ((*rtree)->root == NULL) return;

   rtree_free_node(*rtree, (*rtree)->root, 0, (*rtree)->depth); 
   (*rtree)->root = 0;   /* rtree_free_node already frees mem */

   rtree_mem_free_tree (rtree); 
}

int rtree_empty (rtree_const_handle tree)
{

    /* null tree must be empty, right?*/
    if (!tree) return 1;
    return (rtree_get_count (tree) == 0);
}

void rtree_clear (rtree_handle rtree)
{
   if (!rtree->count)
      return; 

   assert (rtree->root); 

   rtree_free_node (rtree, rtree->root, 0, rtree->depth); 

   /* rtree_free_node also free'd the root node itself,
    * and we want every rtree (also an empty on) to have a valid rtree node
    * so create a new one */

   rtree->root = rtree_mem_alloc_node();
   rtree->root->range.start = 0; 
   rtree->root->range.stop = 0; 

   rtree->count = 0; 
   rtree->depth = 0; 
}


rtree_handle rtree_create ()
{
   int i; 
   rtree_handle ret = rtree_mem_alloc_tree(); 
   assert (ret); 

   ret->depth = 0; 
   ret->root = rtree_mem_alloc_node();
   ret->root->range.start = 0; 
   ret->root->range.stop = 0; 

   ret->count = 0; 

   for (i=0; i<RTREE_CHILD_MAX; ++i)
      ret->root->entry[i]=0; 

   rtree_set_splitfunc (ret, 0, 0); 
   rtree_set_freefunc (ret, 0, 0); 

   return ret; 
}

/*
 * Find the (possibly leaf) node at depth wanted that will cause the least extension 
 * of the bounding box when adding the specified range
 */
static rtree_node * rtree_add_choosenode (rtree_node * node, const rtree_range * range,
      int depth, int treedepth, int wanted)
{
   int i;
   int childcount = RTREE_CHILD_MAX; 
   int addchild; 

   RTREE_RANGE_TYPE addrange; 
   RTREE_RANGE_TYPE addincrease; 
   RTREE_RANGE_TYPE tmprange; 
   RTREE_RANGE_TYPE increase[RTREE_CHILD_MAX];

   /* cannot be looking below desired depth */ 
   assert (depth <= wanted); 

   if (depth == wanted)
      return node; 

   /* calculate increase of the node when adding the rectangle there */
   for (i=0; i<RTREE_CHILD_MAX; ++i)
   {
      if (!node->child[i])
      {
	 childcount = i; 
	 break;
      }
      increase[i] = rtree_range_calc_extension (&node->child[i]->range, range);
   }

   /* select the one with the least increase */
   addchild = 0; 
   addincrease = increase[0]; 
   addrange = rtree_range_size(&node->child[addchild]->range); 

   assert (childcount > 0); 
   for (i=1; i<childcount; ++i)
   {
      if (increase[i] > addincrease)
	 continue;

      /* if the increase is the same, keep the one that
       * has the smallest range */
      if (increase[i] == addincrease)
      {
	 tmprange = rtree_range_size(&node->child[i]->range); 
	 if (tmprange >= addrange)
	    continue; 
      }

      addchild=i; 
      addrange =tmprange;
      addincrease=increase[addchild]; 
   }

   /* we found the best child node to insert the new range */
   return rtree_add_choosenode (node->child[addchild], range, 
	 depth+1, treedepth, wanted); 
}

/* Partition the ranges from source in two sections; 
 * Denoted by 0 or 1 in mapping
 * source and mapping are of size count;
 * All entries need to be checked [there could be 0 pointers in the middle]
 */
static void rtree_add_splitnode_decide (const rtree_range ** source,
      int * mapping, int count)
{
   int i; 
   RTREE_RANGE_TYPE mostleft; 
   RTREE_RANGE_TYPE mostright; 
   int mostrightid, mostleftid; 

   assert (source); 
   assert (source[0]); 

   mostleft = source[0]->start; 
   mostright = source[0]->stop; 
   mostleftid = mostrightid = 0; 

   /* find the points most to the left and to the right */
   for (i=0; i<count; ++i)
   {
      if (!source[i])
	 continue; 
      
      if (source[i]->start < mostleft)
      {
      	mostleft = source[i]->start; 
	mostleftid = i; 
      }
     
      if (source[i]->stop > mostright)
      {
	 mostright = source[i]->stop; 
	 mostrightid = i; 
      }
   }

   if (mostleftid == mostrightid)
   {
      fprintf(stderr, "rtree: complete overlap detected; node splitting"
	    " doesn't handle this well!\n"); 
      /* select another one */
      if (mostleftid == (count-1) || !source[mostleftid+1])
	 mostrightid=0; 
      		else
	 mostrightid=mostleftid+1; 
   }

   assert (mostleftid != mostrightid); 
   assert (source[mostleftid]);
   assert (source[mostrightid]); 
   mapping[mostleftid]=0;
   mapping[mostrightid]=1; 

   /* go over rest and decide which side they belong on */
   for (i=0; i<count; ++i)
   {
      RTREE_RANGE_TYPE leftdiff;
      RTREE_RANGE_TYPE rightdiff; 

      if (!source[i])
	 mapping[i]=0;

      if (i == mostleftid || i == mostrightid)
	 continue; 

      leftdiff = source[i]->stop - source[mostleftid]->stop; 
      rightdiff = source[mostrightid]->start - source[i]->start;  
      mapping[i] = (leftdiff > rightdiff ? 1 : 0); 
   }
}


/*
 * Collect range pointers from the node
 */
static void rtree_add_getranges (const rtree_node * source, const rtree_range ** ranges, int leaf)
{
   int i=0; 
   if (leaf)
   {
      for (i=0; i<RTREE_CHILD_MAX; ++i)
	 *ranges++ = &source->entry[i]->range; 
   }
   else
   {
      for (i=0; i<RTREE_CHILD_MAX; ++i)
	 *ranges++ = &source->child[i]->range; 
   }
}


/* 
 * Fix parent pointer for non-leaf nodes
 */
static void rtree_node_fixparent (rtree_node * node)
{
   int i;
   for (i=0; i<RTREE_CHILD_MAX; ++i)
   {
      if (!node->child[i])
	 break;
      node->child[i]->parent = node; 
   }
}

/* split node [nonleaf]
 * Does not update bounding range
 * Sets parent pointer and the parent pointer of the children
 */
static void rtree_add_splitnode_node (rtree_const_handle tree, 
      rtree_node ** source, rtree_node ** s1, 
      rtree_node ** s2, rtree_node * newchild)
{
   int i; 
   int mapping[RTREE_CHILD_MAX+1]; 
   const rtree_range * ranges[RTREE_CHILD_MAX+1]; 

   int leftpos  = 0;
   int rightpos = 0; 

   assert (newchild);

   ranges[RTREE_CHILD_MAX] = &newchild->range; 
   rtree_add_getranges (*source, ranges,FALSE);   

   tree->splitfunc (ranges, mapping, RTREE_CHILD_MAX+1, 
	 tree->splitdata); 

   *s1 = rtree_mem_alloc_node (); 
   *s2 = rtree_mem_alloc_node (); 

   (*s1)->parent = (*source)->parent;
   (*s2)->parent = (*source)->parent; 

   /* add old children */ 
   for (i=0; i<RTREE_CHILD_MAX; ++i)
   {
      if (mapping[i])
	 (*s2)->child[rightpos++] = (*source)->child[i]; 
   		else
	 (*s1)->child[leftpos++] = (*source)->child[i];

      assert (mapping[i] == 0 || mapping[i] ==1 ); 
   }

   /* add new child */
   if (mapping[RTREE_CHILD_MAX])
   {
      
      (*s2)->child[rightpos++] = newchild; 
#ifdef RTREE_SORT_ENTRIES
      rtree_node_sort (*s2, rightpos, 0); 
#endif
   }
   	else
   {
      (*s1)->child[leftpos++] = newchild; 
#ifdef RTREE_SORT_ENTRIES
      rtree_node_sort (*s1, leftpos, 0); 
#endif
   }

   /* Fix parent pointer in child nodes */ 
   rtree_node_fixparent (*s1);
   rtree_node_fixparent (*s2); 

   rtree_mem_free_node (source);  
}

/* split node [must be leaf] 
 * Does not update bounding range
 * Sets parent pointer to that of the source node
 **/
void rtree_add_splitnode_leaf (rtree_const_handle tree, 
      rtree_node ** source, rtree_node ** s1, rtree_node ** s2, 
      rtree_entry * newentry)
{
   int i; 
   int mapping[RTREE_CHILD_MAX+1]; 
   const rtree_range * ranges[RTREE_CHILD_MAX+1]; 

   int leftpos  = 0;
   int rightpos = 0; 

   assert (newentry);

   ranges[RTREE_CHILD_MAX] = &newentry->range; 
   rtree_add_getranges (*source, ranges,TRUE);   

   /* calculate new distribution */
   tree->splitfunc (ranges, mapping, RTREE_CHILD_MAX+1, tree->splitdata); 

   *s1 = rtree_mem_alloc_node (); 
   *s2 = rtree_mem_alloc_node (); 

   (*s1)->parent = (*source)->parent;
   (*s2)->parent = (*source)->parent; 

   /* add old entries */ 
   for (i=0; i<RTREE_CHILD_MAX; ++i)
   {
      if (mapping[i])
	 (*s2)->entry[rightpos++] = (*source)->entry[i]; 
   		else
	 (*s1)->entry[leftpos++] = (*source)->entry[i];

      assert (mapping[i] == 0 || mapping[i] ==1 ); 
   }

   /* add new entry */
   if (mapping[RTREE_CHILD_MAX])
   {
      (*s2)->entry[rightpos++] = newentry; 
#ifdef RTREE_SORT_ENTRIES
      rtree_node_sort (*s2, rightpos, 1); 
#endif
   }
   	else
   {
      (*s1)->entry[leftpos++] = newentry; 
#ifdef RTREE_SORT_ENTRIES
   rtree_node_sort (*s1, leftpos, 1); 
#endif
   }



   /* entry nodes have no parent pointer so nothing to correct */ 
   rtree_mem_free_node (source);  
}


static void rtree_node_fix_extent (rtree_node * node, int leaf)
{
   int i=0; 
   assert (node); 

   /* could do some define hack to avoid duplicate code */

   if (!leaf)
   {
      if (node->child[0])
      {
	 node->range = node->child[0]->range; 
      }
      else
      {
	 node->range.start = node->range.stop = 0; 
      }

      for (i=1; i<RTREE_CHILD_MAX; ++i)
      {
	 if (!node->child[i])
	    break;
	 rtree_range_extend (&node->range, &node->child[i]->range); 
      }
   }
   else
   {
      if (node->entry[0])
      {
	 node->range = node->entry[0]->range; 
      }
      else
      {
	 node->range.start = node->range.stop = 0; 
      }

      for (i=1; i<RTREE_CHILD_MAX; ++i)
      {
	 if (!node->entry[i])
	    break;
	 rtree_range_extend (&node->range, &node->entry[i]->range); 
      }
   }
}

/* find an empty slot in the node; RTREE_CHILD_MAX if full */ 
static inline int rtree_node_findempty_child (const rtree_node * n)
{
   int i;
   assert (n); 
   for (i=0; i<RTREE_CHILD_MAX; ++i)
   {
      if (!n->child[i])
	 break; 
   }
   return i; 
}

/* find empty entry */
static inline int rtree_node_findempty_entry (const rtree_node * n)
{
   int i;
   assert (n); 
   for (i=0; i<RTREE_CHILD_MAX; ++i)
   {
      if (!n->entry[i])
	 break; 
   }
   return i; 
}

/*
 * Fixes extent of n1 and n2 if not zero; 
 * If n2 is non-zero, add it to the parent of n1 (which could result in node
 * splits)
 */
static void rtree_add_adjusttree (rtree_handle tree, rtree_node * n1, 
      rtree_node * n2, int depth)
{
   rtree_node * s1 = 0 ;
   rtree_node * s2 = 0; 

   assert (n1); 
   assert (depth >= 0); 

   /* adjust bounding box */
   rtree_node_fix_extent (n1, depth==tree->depth);

   if (n2)
      rtree_node_fix_extent (n2, depth==tree->depth); 

   /* Is n1 the root? */
   if (!n1->parent)
   {
      assert (!n2 || !n2->parent);

      if (n2)
      {
	 /* split went up to the root, create new root node */
	 rtree_node * newroot = rtree_mem_alloc_node ();

	 newroot->parent = 0;
	 newroot->child[0]=n1; 
	 newroot->child[1]=n2; 

	 /* fix parent pointer */ 
	 rtree_node_fixparent (newroot); 
	 rtree_node_fix_extent (newroot, 0 == tree->depth); 

	 tree->root = newroot; 
	 ++(tree->depth);       
      }

      /* extent of n1, n2 and parent is ok; we're done */
      return; 
   }

   /* was there a split? */
   if (n2)
   {
      /* try to add new node to the parent of n1 */
      rtree_node * parent;
      int pos; 
      
      /* try to add the new node to the parent  */ 
      /* the parent is for sure a normal non-leaf node */
      parent = n1->parent; 
      assert (n2->parent == parent); 
      pos = rtree_node_findempty_child (parent); 
      
      if (pos == RTREE_CHILD_MAX)
      {
	 const rtree_node * old = parent;
	 rtree_node * parentparent = parent->parent; 
	 int j; 

	 /* Parent node is also full; split */
	 rtree_add_splitnode_node (tree, &parent, &s1, &s2, n2);  

	 if (parentparent)
	 {
	    /* parent gets freed because of splitnode; but the parent of
	     * parent has a child pointer to parent, which needs to be updated
	     */
	    for (j=0; j<RTREE_CHILD_MAX; ++j)
	    {
	       if (parentparent->child[j]==old)
	       {
		  parentparent->child[j] = s1; 
		  break; 
	       }
	    }
	    assert (j!=RTREE_CHILD_MAX); 
	 }
      }
      else
      {
	 /* room in node; add */
	 assert (!parent->child[pos]); 
	 assert (n2->parent == parent); 
	 parent->child[pos] = n2; 
	 
	 s1 = n1->parent; 
	 s2 = 0; 
#ifdef RTREE_SORT_NODES
         rtree_node_sort (parent, pos+1, 0); 
#endif
      }
   }
   else
   {
      s1 = n1->parent; 
      s2 = 0; 
   }
   
   /* adjust the parent */
   rtree_add_adjusttree (tree, s1, s2, depth-1);  
}

/* try to add entry to the given node 
 * Does NOT update bounding rectangle 
 * Return true if succeeded
 * */
static int rtree_add_try_leaf (rtree_node * node, rtree_entry * new) 
{
   int i; 
   for (i=0; i<RTREE_CHILD_MAX; ++i)
   {
      if (!node->entry[i])
      {
	 /* found empty space */
	 break; 
      }
   }

   if (i != RTREE_CHILD_MAX)
   {
	 node->entry[i] = new; 
#ifdef RTREE_SORT_ENTRIES
         /* for now, just resort the whole thing afterwards */
         rtree_node_sort (node, i+1, 1); 
#endif
         return 1; 
   }
   else
      return 0; 
}

/* Add entry; Does NOT update tree->count  */ 
static void rtree_add_entry (rtree_handle tree, 
      rtree_entry * newentry)
{
   rtree_node * addpoint = 0; 
   rtree_node * s1 = 0; 
   rtree_node * s2 = 0; 
   int done; 

   /* find a leaf node where the new range fits best */
   addpoint = rtree_add_choosenode (tree->root, &newentry->range, 
	 0, tree->depth, tree->depth); 

   /* if there is still space in the node, add it there
    * This does NOT update the bounding box of the node */
   done = rtree_add_try_leaf (addpoint, newentry); 

   if (done)
   {
      /* add worked */
      s1 = addpoint; s2 = 0; addpoint = 0; 
   }
   else
   {
      int j; 
      rtree_node * parentparent = addpoint->parent; 
      rtree_node * old = addpoint; 

      /* leaf node is full; Split and add new value to one of the nodes */ 
      /* we know this is a leaf node */
      rtree_add_splitnode_leaf (tree, &addpoint, &s1, &s2, newentry); 

      /* if the node has a parent, we need to update the child pointer since
       * 'addpoint' node could be freed during the split */
      if (parentparent)
      {
	 for (j=0; j<RTREE_CHILD_MAX; ++j)
	 {
	    if (parentparent->child[j]==old)
	    {
	       parentparent->child[j] = s1; 
	       break; 
	    }
	 }
	 assert (j!=RTREE_CHILD_MAX); 
      }
   }
   
   /* now addpoint is 0
    * It there was a split, s2 is not null; Together s1 and s2 contain 
    * all entries of addpoint and the newly added entry
    * Adjust bounding boxes along the path to the root; 
    * Integrate the splitted node if needed
    */
   rtree_add_adjusttree (tree, s1, s2, tree->depth);

}

/* try to add newnode to the given node 
 * Does NOT update bounding rectangle 
 * Return true if succeeded
 * */
static int rtree_add_try_node (rtree_node * node, rtree_node * new) 
{
   int i; 
   for (i=0; i<RTREE_CHILD_MAX; ++i)
   {
      if (!node->child[i])
      {
	 /* found empty space */
	 break; 
      }
   }

#ifndef RTREE_NODE_SORT
   if (i != RTREE_CHILD_MAX)
   {
	 node->child[i] = new; 
         return 1; 
   }
   else
      return 0; 
#else
   /* there is space, find the best position and add there */
   assert (false); 
#endif
}


/* Try to add newnode to node 'node'  */ 
static void rtree_add_node (rtree_handle tree, rtree_node * node, rtree_node * newnode)
{
   rtree_node * s1; 
   rtree_node * s2; 
   int done; 

   /* if there is still space in the node, add it there
    * This does NOT update the bounding box of the node */
   done = rtree_add_try_node (node,newnode); 

   if (done)
   {
      /* add worked */
      s1 = node; s2 = 0; node = 0; 
   }
   else
   {
      int j; 
      rtree_node * parentparent = node->parent; 
      rtree_node * old = node; 

      /* leaf node is full; Split and add new value to one of the nodes */ 
      /* we know this is a leaf node */
      rtree_add_splitnode_node (tree, &node, &s1, &s2, newnode); 

      /* if the node has a parent, we need to update the child pointer since
       * 'addpoint' node could be freed during the split */
      if (parentparent)
      {
	 for (j=0; j<RTREE_CHILD_MAX; ++j)
	 {
	    if (parentparent->child[j]==old)
	    {
	       parentparent->child[j] = s1; 
	       break; 
	    }
	 }
	 assert (j!=RTREE_CHILD_MAX); 
      }
   }
   
   /* now addpoint is 0
    * It there was a split, s2 is not null; Together s1 and s2 contain 
    * all entries of addpoint and the newly added entry
    * Adjust bounding boxes along the path to the root; 
    * Integrate the splitted node if needed
    */
   rtree_add_adjusttree (tree, s1, s2, tree->depth);

}


void rtree_add (rtree_handle tree, const rtree_range * range, 
      RTREE_DATA_TYPE data)
{
   rtree_entry * newentry = 0; 

   /* Create entry for it */
   newentry = rtree_mem_alloc_entry (range, data); 

   /* add it */ 
   rtree_add_entry (tree, newentry); 
   
   /* update count */ 
   ++tree->count; 
}




int rtree_walk_internal (const rtree_node * node, rtree_callback callback, 
      void * extra, int depth, int treedepth )
{
   int i; 

   assert (node); 
   if (depth < treedepth)
   {
      for (i=0; i<RTREE_CHILD_MAX; ++i)
      {
	 if (!node->child[i])
	 {
	    assert (i>0); 
	    break; 
	 }
	 if (!rtree_walk_internal (node->child[i], 
		  callback, extra, depth+1, treedepth))
	    return FALSE; 
      }
      return TRUE; 
   }

   /* lowest level */ 
   for (i=0; i<RTREE_CHILD_MAX; ++i)
   {
      if (!node->entry[i])
	break; 

      if (!callback (&node->entry[i]->range, 
	       &node->entry[i]->data, extra))
	 return FALSE;
   }
   return TRUE; 
}

int rtree_walk (rtree_const_handle tree, rtree_callback callback, void * extra)
{
    if(tree == NULL) return 1;
    if (tree->root == NULL) return 1;

   return rtree_walk_internal (tree->root, callback, extra, 0, tree->depth); 
}

static int rtree_walk_all_internal (const rtree_node * node, rtree_callback_all callback, 
      rtree_callback_all_info * info, int depth, int treedepth )
{
   int i; 

   assert (node); 
	

   /* show the current node */ 
   info->depth = depth; 
   assert (info->treedepth == treedepth); 
   info->nodeid = (void *) node;
   info->parentid = node->parent;
   info->data = 0; 
   info->range = &node->range; 

   if (!callback (info))
	    return FALSE; 

   if (depth < treedepth)
   {
      for (i=0; i<RTREE_CHILD_MAX; ++i)
      {
	 if (!node->child[i])
	 {
	    assert (i>0); 
	    break; 
	 }

	 if (!rtree_walk_all_internal (node->child[i], 
		  callback, info, depth+1, treedepth))
	    return FALSE; 
      }
      return TRUE; 
   }

   /* lowest level */
   info->parentid = (void *) node; 
   info->depth = depth + 1; 
   for (i=0; i<RTREE_CHILD_MAX; ++i)
   {
      if (!node->entry[i])
	break; 

      /* show the entries of the leaf node; depth will be > treedepth */ 
      info->nodeid = node->entry[i]; 
      info->range = &node->entry[i]->range; 
      info->data = &node->entry[i]->data; 
      if (!callback (info))
	 return FALSE;
   }
   return TRUE; 
}


int rtree_walk_all (rtree_const_handle tree, rtree_callback_all callback, void * extra)
{
   assert (tree); 
   assert (tree->root); 

   rtree_callback_all_info info; 
   info.extra =extra; 
   info.treedepth = tree->depth; 
   info.tree = tree; 


   return rtree_walk_all_internal (tree->root, callback, &info, 0, tree->depth); 
}


void rtree_get_range (rtree_const_handle tree, rtree_range  * range)
{
   assert (range); 
   assert (tree);
   assert (tree->root); 
   *range = tree->root->range; 
}

int rtree_get_depth (rtree_const_handle tree)
{
   assert (tree); 
   return tree->depth; 
}


static inline void rtree_dump_internal_indent (int amount)
{
   int i;
   for (i=0; i<amount; ++i)
      printf("  "); 
}


static int rtree_dumpfunc (rtree_const_handle tree, const rtree_node * node, int depth, 
      void * extra)
{
   int i;

   rtree_dump_internal_indent (2*depth); 
   printf ("[" RTREE_RANGE_TYPE_PRINTF "," RTREE_RANGE_TYPE_PRINTF "[\n", 
         node->range.start, node->range.stop); 

   if (depth == tree->depth)
   {
      for (i=0; i<RTREE_CHILD_MAX; ++i)
      {
	 if (!node->entry[i])
	    break; 
	 rtree_dump_internal_indent (2*depth+1);
	 printf ("=> ENTRY [" RTREE_RANGE_TYPE_PRINTF "," RTREE_RANGE_TYPE_PRINTF "[ "
               RTREE_DATA_TYPE_PRINTF "\n", node->entry[i]->range.start, 
	       node->entry[i]->range.stop, node->entry[i]->data); 
      }
   }
   return TRUE; 
}

static int rtree_visit_nodes_internal (rtree_const_handle tree, const rtree_node * node, 
      int depth, rtree_callback_visit func, void * extra)
{
   int i; 
   if (!func (tree, node, depth, extra))
      return FALSE; 

   if (depth < tree->depth)
   {
      for (i=0; i<RTREE_CHILD_MAX; ++i)
      {
	 if (!node->child[i])
	    break; 
         if (!rtree_visit_nodes_internal (tree, node->child[i], depth+1, 
                  func, extra))
            return FALSE; 
      }
   }
   return TRUE; 
}

/*
 * For every node, call func with extra
 * Stops at the leaf nodes 
 */
static int rtree_visit_nodes (rtree_const_handle tree, rtree_callback_visit func, void * extra)
{
   return rtree_visit_nodes_internal (tree, tree->root, 0, func, extra);  
}

void rtree_dump (rtree_const_handle tree)
{
   assert (tree);
   rtree_visit_nodes (tree, rtree_dumpfunc, 0); 
}

int rtree_get_child_min (rtree_const_handle tree)
{
   assert (tree); 
   return RTREE_CHILD_MIN; 
}

int rtree_get_child_max (rtree_const_handle tree)
{
   assert (tree); 
   return RTREE_CHILD_MAX; 
}

int rtree_get_count (rtree_const_handle tree)
{
   assert (tree); 
   return tree->count; 
}

void rtree_set_freefunc (rtree_handle tree, rtree_callback func, 
      void * extra)
{
   assert (tree); 
   tree->freedata = extra;
   tree->freefunc = func; 
}


static void rtree_default_splitfunc (const rtree_range ** sources, 
      int * mapping, int count, void * extra)
{
    rtree_add_splitnode_decide(sources,mapping,count); 
}

void rtree_set_splitfunc (rtree_handle tree, rtree_callback_split func, 
      void * extra)
{
   if (!func)
   {
      tree->splitdata = 0; 
      tree->splitfunc = rtree_default_splitfunc; 
      return; 
   }
   tree->splitdata = extra; 
   tree->splitfunc = func; 
}

static int rtree_validate_node (rtree_const_handle tree, 
      const rtree_node * node, int depth, void * extra)
{
   const rtree_range * ranges[RTREE_CHILD_MAX]; 
   int i; 
   rtree_range check; 

   assert (node);

   if (depth == tree->depth)
   {
      for (i=0; i<RTREE_CHILD_MAX; ++i)
         ranges[i] = (node->entry[i] ? &node->entry[i]->range : 0); 
   }
   else
   {
      for (i=0; i<RTREE_CHILD_MAX; ++i)
         ranges[i] = (node->child[i] ? &node->child[i]->range : 0); 
   }

   /* Check that node->range equals the extent of 
    * the ranges of the children */
   if (ranges[0])
   {
      check = *ranges[0]; 
      for (i=1; i<RTREE_CHILD_MAX; ++i)
      {
	 if (!ranges[i])
	    break;
	 rtree_range_extend (&check, ranges[i]); 
      }
   }
   else
   {
      fprintf (stderr, "RTree error: empty node!\n"); 
      return FALSE; 
   }

   if (rtree_range_equals (&check, &node->range))
      return TRUE; 
  
   fprintf (stderr, "Error in node range!\n"); 
   assert (FALSE); 
   return FALSE; 
}

static int rtree_node_count_children (const rtree_node * node, int leaf)
{
   int i=0; 
   if (leaf)
   {
      for (i=0; i<RTREE_CHILD_MAX; ++i)
         if (!node->entry[i])
            break;
   }
   else
   {
      for (i=0; i<RTREE_CHILD_MAX; ++i)
         if (!node->child[i])
            break; 
   }
   return i; 
}


/*
 * Returns RTREE_CHILD_MAX if the child could not be found in the node */ 
static inline int rtree_node_find_child (const rtree_node * n, const rtree_node * child)
{
   int i;
   assert (n); 
   for (i=0; i<RTREE_CHILD_MAX; ++i)
   {
      if (n->child[i]==child)
	 break; 
   }
   return i; 
}


static int rtree_validate_count (rtree_const_handle tree, const rtree_node * node, 
      int depth, void * extra)
{
   int * count = (int *) extra; 

   assert (node); assert (extra); 
   
   /* only interested in leaf nodes */ 
   if (depth != tree->depth)
      return TRUE; 

   *count += rtree_node_count_children (node, TRUE); 
   return TRUE; 
}

static int rtree_validate_parent (rtree_const_handle tree, const rtree_node * node, 
      int depth, void * extra)
{
   int pos; 

   assert (node); 

   if (!depth)
   {
      /* root node */ 
      assert (!node->parent); 
      return (!node->parent); 
   }

   pos = rtree_node_find_child (node->parent, node); 
   if (pos == RTREE_CHILD_MAX)
   {
      fprintf (stderr, "RTree: node->parent link incorrect!\n");
      assert (node->parent->child[pos] == node); 
      return FALSE; 
   }

   return TRUE; 
}

/*
 * Validate tree structure
 *    - test node extent
 */
int rtree_check (rtree_const_handle tree)
{
   int count = 0; 

   /* check extent for the nodes */ 
   if (!rtree_visit_nodes (tree, rtree_validate_node, 0))
         return FALSE;

   /* check child count consistency */
   rtree_visit_nodes (tree, rtree_validate_count, (void*)&count);
   assert (count == tree->count); 
   if (count != tree->count)
   {
      fprintf (stderr, "RTree: tree->count not consistent with "
            "entry count\n"); 
      return FALSE; 
   }

   /* check parent pointer */ 
   if (!rtree_visit_nodes (tree, rtree_validate_parent, 0))
         return FALSE;

   return TRUE; 
}


/*************** node removal ******************/

static void rtree_node_remove_entry (rtree_node * node, int entry)
{
   int pos = entry; 

   assert (node); 
   assert (node->entry[entry]); 
   
   node->entry[entry] = 0; 


#ifdef RTREE_LEAF_SORT
   /***** INCORRECT !!  [ last pointer twice in index !  ] *****/ 
   /* move up the remaining entries */
   for (pos = entry; pos<RTREE_CHILD_MAX-1; ++pos)
   {
      node->entry[pos] = node->entry[pos+1]; 
      if (!node->entry[pos])
         break; 
   }
   assert (false); 
#else
   /* Find first zero pointer after the removed item */
   for (pos=entry+1; pos<RTREE_CHILD_MAX; ++pos)
   {
      if (!node->entry[pos])
         break; 
   }
   /* now pos = last_not_null_idx + 1; If equal to pos+1 -> no elements left */
   assert ((pos <= RTREE_CHILD_MAX) && pos > 0);
   if (pos != entry+1)
   {
      node->entry[entry] = node->entry[pos-1]; 
      node->entry[pos-1] = 0; 
   }
#endif

   rtree_node_fix_extent (node, TRUE); 
}


static void rtree_node_remove_child (rtree_node * node, int childnum)
{
   int pos = childnum; 

   assert (node); 
   assert (node->child[childnum]); 
   
   node->child[childnum] = 0; 

#ifdef RTREE_NODE_SORT
   /* move up the remaining entries */
   for (pos = childnum; pos<RTREE_CHILD_MAX-1; ++pos)
   {
      node->child[pos] = node->child[pos+1]; 
      if (!node->child[pos])
         break; 
   }
   assert (false); 
#else
   /* put last child in removed slot*/
   for (pos=childnum+1; pos<RTREE_CHILD_MAX; ++pos)
   {
      if (!node->child[pos])
         break; 
   }
   assert ((pos <= RTREE_CHILD_MAX) && pos > 0);
   if (pos != childnum+1)
   {
      node->entry[childnum] = node->entry[pos-1]; 
      node->entry[pos-1] = 0; 
   }
#endif

   rtree_node_fix_extent (node, FALSE); 
}

static void rtree_condensetree (rtree_handle tree, rtree_node * node, int depth)
{
   int remove; 
   rtree_node * parent = 0;
   int childcount; 

   assert (node);

   parent = node->parent;

   /* check if node needs to be removed */
   childcount = rtree_node_count_children(node, depth==tree->depth); 

   remove =(childcount < rtree_get_child_min(tree) ? 1 : 0)
                && depth; /* cannot remove root node */

   if (remove)
   {
      int childnum; /* child number of this child in the parent */

      /* not enough children in the node, remove the node 
       * which means also removing something from the parent node
       */

      if (parent)
      {
         /* possible performance improvement: avoid the search in the parent
          * node to find our pointer */

         childnum = rtree_node_find_child(parent, node);
         assert (childnum != RTREE_CHILD_MAX); 

         rtree_node_remove_child (parent, childnum); 
      }
   }
        else
   {
      /* the node stays, just fix the extent */ 
      rtree_node_fix_extent (node, depth == tree->depth); 
   }
       
   /* if we remove or not, need to go up to the parent */
   if (parent)
   {
      assert (depth); 
      rtree_condensetree (tree, node->parent, depth-1); 
   }

   if (remove)
   {
      /* if we removed a node, insert its remaining children again at the same level */
      /* then free the node */ 
      int i; 

      for (i=0; i<RTREE_CHILD_MAX; ++i)
      {
         if (depth == tree->depth)
         {
            if (!node->entry[i])
               break; 

            /* reinserting leaf node */ 
            rtree_add_entry (tree, node->entry[i]); 
         }
         else
         {
            rtree_node * insertpoint; 
            /* Reinserting a non-leaf node */
            if (!node->child[i])
               break; 

            /* depth is the depth of this node */ 
            insertpoint = rtree_add_choosenode (tree->root, &node->child[i]->range, 
	                0, tree->depth, depth); 
               
            rtree_add_node (tree, insertpoint, node->child[i]); 
         }
      }
      rtree_mem_free_node (&node); 
   }
}

/* See if we can reduce the height of the tree */ 
static void rtree_remove_checkroot (rtree_handle tree)
{
   rtree_node * tmp;
   int childcount;


   /* the root node has is a leaf node if the depth of the tree is 0 */ 
   childcount = rtree_node_count_children (tree->root, tree->depth == 0); 
   assert (childcount || !tree->depth); 

   /* if the there is only one node we cannot remove is (depth=0) */
   if (childcount > 1 || !tree->depth)
      return; 

   tmp = tree->root->child[0]; 

   assert (tmp); 
   assert (rtree_range_equals (&tree->root->range, 
	    &tmp->range)); 

   rtree_mem_free_node (&tree->root); 
   tree->root = tmp; 
   assert (tree->depth); 
   --tree->depth; 
   tree->root->parent = 0; 
}

int rtree_remove (rtree_handle tree, const rtree_range * range, 
      RTREE_DATA_TYPE * data)
{
   rtree_node * nodeptr = 0;
   rtree_entry * entry = 0; 
   int entrynum = -1; 

   /* try to find node first */ 
   rtree_find_internal (tree->root, range, 0, tree->depth, 
         (const rtree_node **) /* const cast */ &nodeptr, &entrynum);

   if (!nodeptr)
      return FALSE;

   /* found, store data */
   entry = nodeptr->entry[entrynum]; 
   if (data)
      *data = entry->data; 

   /* call free function if there is one */
   if (tree->freefunc)
      tree->freefunc (&entry->range, &entry->data, tree->freedata);

   /* remove entry pointer from node */
   rtree_node_remove_entry (nodeptr, entrynum); 

   /* remove entry from tree */ 
   rtree_condensetree (tree, nodeptr, tree->depth); 

   /* if the root has only one child decrease the height of the tree */ 
   rtree_remove_checkroot (tree); 

   /* free entry memory */
   rtree_mem_free_entry (&entry);
   
   /* adjust item count */
   --tree->count; 

   return TRUE; 
}

static rtree_entry * rtree_copy_entry (rtree_entry * entry, 
      rtree_callback_copy copy, void * extra)
{
   rtree_entry * newentry = rtree_mem_alloc_entry (&entry->range,
         0); 

   if (copy)
      copy (&entry->range, &entry->data, &newentry->data);
        else
      newentry->data = entry->data; 

   return newentry; 
}

static rtree_node * rtree_copy_internal (const rtree_node * node, 
      int depth, int maxdepth, rtree_callback_copy copy, void * extra)
{
   int i; 
   rtree_node * newnode; 
   assert (node); 

   newnode = rtree_mem_alloc_node (); 
   assert (newnode); 

   newnode->range = node->range; 
   newnode->parent = 0; 

   if (depth == maxdepth)
   {
      /* copy entries */ 
      for (i=0; i<RTREE_CHILD_MAX; ++i)
      {
         if (node->entry[i])
            newnode->entry[i] = rtree_copy_entry (node->entry[i], 
                  copy, extra); 
          else
            newnode->entry[i] = 0; 
         /* entries have no parent pointer */ 
      }
   }
   else
   {
      /* copy children */ 
      for (i=0; i<RTREE_CHILD_MAX; ++i)
      {
         if (node->child[i])
         {
            newnode->child[i] = rtree_copy_internal (node->child[i], 
               depth+1, maxdepth, copy, extra);
            /* fix parent pointer */
            newnode->child[i]->parent = newnode; 
         }
               else
            newnode->child[i] = 0;  
      }
   }
   return newnode; 
}

rtree_handle rtree_copy (rtree_const_handle tree, rtree_callback_copy copy, 
      void * extra)
{
   rtree_handle newtree = rtree_mem_alloc_tree (); 
   *newtree = *tree; 

   newtree->root = rtree_copy_internal (tree->root, 0, tree->depth, copy, extra); 
   assert (0 == tree->root->parent); 
   return newtree; 
}

/*=========================================================================
 * RTree range helper functions                                                       *
 *=========================================================================*/
static inline int rtree_range_empty (const rtree_range * r)
{
   return (r->start >= r->stop); 
}

static inline RTREE_RANGE_TYPE rtree_range_size (const rtree_range * r)
{
   assert (r->start <= r->stop); 
   return (r->stop - r->start); 
}

/* calculate how much r1 would extend when adding r2 to it */
static inline RTREE_RANGE_TYPE rtree_range_calc_extension (const rtree_range * r1, 
      const rtree_range * r2)
{
   rtree_range d; 
   
   rtree_range_extent (r1, r2, &d); 
   return rtree_range_size(&d) - rtree_range_size(r1); 
}

static inline RTREE_RANGE_TYPE rtree_range_type_min 
		(RTREE_RANGE_TYPE a, RTREE_RANGE_TYPE b)
{
   return (a < b ? a : b); 
}

static inline RTREE_RANGE_TYPE rtree_range_type_max 
	(RTREE_RANGE_TYPE a, RTREE_RANGE_TYPE b)
{
   return (a > b ? a : b); 
}


/* Switch pointers if needed so that range1 comes first */
static inline void rtree_range_order (const rtree_range * * range1, 
      const rtree_range * * range2)
{
   if ((*range1)->start > (*range2)->start)
   {
      const rtree_range * tmp = *range1; 
      *range1 = *range2; *range2 = tmp; 
   }
}

static inline int rtree_range_has_overlap (const rtree_range * range1, 
      		const rtree_range * range2)
{
   return (rtree_range_type_max(range1->start, range2->start)
	 	< rtree_range_type_min(range1->stop, range2->stop)); 
}

static inline int rtree_range_shared (const rtree_range * range1, 
      const rtree_range * range2, rtree_range * dest)
{
   dest->start = rtree_range_type_max(range1->start, range2->start); 
   dest->stop = rtree_range_type_min(range1->stop, range2->stop);
   return (dest->start < dest->stop); 
}

static inline void rtree_range_extent (const rtree_range * range1, 
      const rtree_range * range2, rtree_range * dest)
{
   dest->start = rtree_range_type_min(range1->start,range2->start); 
   dest->stop = rtree_range_type_max(range1->stop, range2->stop); 
}

static inline void rtree_range_extend (rtree_range * range1, 
      const rtree_range * range2)
{
   range1->start = rtree_range_type_min(range1->start, range2->start); 
   range1->stop = rtree_range_type_max(range1->stop, range2->stop); 
}


/* check if the ranges are equal */
static inline int rtree_range_equals (const rtree_range * r1, const rtree_range * r2)
{
   assert (r1->start <= r1->stop); 
   if (r1->start == r1->stop)
   {
      return (r2->start == r2->stop); 
   }

   return ((r1->start == r2->start) && (r1->stop == r2->stop)); 
}

/* check if the first range contains all of the second range */
static inline int rtree_range_contains (const rtree_range * r1, const rtree_range * r2)
{
   /* figure out later if two empty ranges contain eachother... */ 
   assert ((r1->start != r1->stop) || (r2->start != r2->stop)); 
   return ((r1->start <= r2->start) && (r1->stop >= r2->stop)); 
}


/****************************************************************************
 * Iterator functions                                                       *
 ****************************************************************************/
static inline void rtree_iterator_validate (rtree_iterator * iter)
{
   assert (iter); 
   assert (iter->tree); 
   assert (iter->createdepth == iter->tree->depth); 
}

/* fast recreate iterator after tree modification */

void rtree_iterator_init (rtree_handle tree, rtree_iterator * iter)
{
   iter->childnum[0] = 0; 
   iter->depth=0; 
   iter->node=tree->root; 
   iter->tree = tree; 
   rtree_iterator_validate (iter); 
}

void rtree_iterator_update (rtree_iterator * iter)
{
   assert (iter); 
   if (iter->createdepth != iter->tree->depth)
   {
      /* realloc faster? */
      free(iter->childnum);
      iter->createdepth = iter->tree->depth; 
      iter->childnum=(int *) malloc (sizeof(int)*iter->createdepth); 
   }

   rtree_iterator_init (iter->tree, iter); 
}

void rtree_iterator_free (rtree_iterator_handle * iter)
{
   assert (iter); 
   free ((*iter)->childnum); 
   rtree_mem_free_iterator (iter); 
}

rtree_iterator_handle rtree_iterator_create (rtree_handle tree)
{
   rtree_iterator * ret = rtree_mem_alloc_iterator (); 
   ret->createdepth = tree->depth; 
   ret->childnum = (int*) malloc(sizeof(int)*ret->createdepth); 
   rtree_iterator_init (tree, ret); 
   rtree_iterator_forward(ret); 
   return ret; 
}

void rtree_iterator_forward (rtree_iterator_handle iter)
{
   /* if we're down and have entries left in the current node,
    * advance child */
   if (iter->depth == iter->tree->depth)
   {
      assert (iter->node); 
      ++iter->childnum[iter->depth]; 
      if (iter->childnum[iter->depth] < 
            rtree_get_child_max(iter->tree) && 
            iter->node->entry[iter->childnum[iter->depth]])
         return; 
   }
   /* go up if needed, then go down the tree */
}

void rtree_iterator_backward (rtree_iterator_handle iter)
{
   assert (FALSE); 
}


