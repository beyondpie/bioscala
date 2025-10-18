package bioscala.LightCoord

import SZUtils.ifelse

import scala.collection.mutable.ListBuffer

/** Genome element range in a given chrom.
  *   - Only start and end position, no chrom info.
  *   - Always from 5' to 3', i.e., startFrom <= endTo
  *   - Follow Bed format, 0-index and range will be
  *     left closed and right open.
  *   - Both start and end positions should be no
  *     smaller than zero.
  *   - No strand information.
  */
type Coord  = (startFrom: Int, endTo: Int)
type Coords = Vector[Coord]

given coordOrd: Ordering[Coord] with {
  def compare(x: Coord, y: Coord): Int = {
    if (x.startFrom < y.startFrom) {
      -1
    } else if (x.startFrom > y.startFrom) {
      1
    } else {
      x.endTo.compare(y.endTo)
    }
  }
}

/** Returns the distance between two coord.
  *   1. Zero if they are overlapped.
  *   2. Positive represents closest gap between them.
  *
  * This funcitons is symmetric.
  *
  * @param x
  *   Coord
  * @param y
  *   Coord
  */
def d(x: Coord, y: Coord): Int = {
  (x.startFrom - y.endTo).max((y.startFrom - x.endTo)).max(0)
}

/** Get the center position of a coord. Here center is
  * startFrom + |__ (endTo - startFrom) / 2 __| The
  * latter part will be integer part of the division.
  * @param x
  *   single coord
  * @return
  */
def getCenter(x: Coord): Coord = {
  val c: Int = x.startFrom + (x.endTo - x.startFrom) / 2
  (startFrom = c, endTo = c + 1)
}

def toStringCoord(g: Coord): String = {
  s"${g.startFrom.toInt}-${g.endTo.toInt}"
}

def mkStringCoord(x: Coord, sep: String = "\t"): String = {
  s"${x.startFrom}${sep}${x.endTo}"
}

/** Check if a vector of genomic coordinates are well
  * sorted or not Ref: StackOverflow about
  * check-whether-a-collection-is-ordered
  * @param x
  * @return
  */
def isCoordSorted(x: Iterable[Coord]): Boolean = {
  x.view.zip(x.tail).forall(((x, y) => coordOrd.lteq(x, y)))
}

/** Intersect between two coords.
  *
  * @param x
  * @param y
  * @return
  */
def intersect(x: Coord, y: Coord): Option[Coord] = {
  if (x.startFrom <= y.startFrom) {
    if (x.endTo <= y.startFrom) {
      None
    } else {
      Option((startFrom = y.startFrom,
              endTo = x.endTo.min(y.endTo)))
    }
  } else {
    if (y.endTo <= x.startFrom) {
      None
    } else {
      Option((startFrom = x.startFrom,
              endTo = y.endTo.min(x.endTo)))
    }
  }
}

/** Check if two Coords overlap (true) or not (false)
  * Oder of input does not matter.
  * @param x
  *   coord
  * @param y
  *   coord
  * @return
  */
def isOvlp(x: Coord, y: Coord): Boolean = {
  ifelse(x.endTo < y.startFrom || x.startFrom >= y.endTo, false,
      true)
}

/** Find overlapped Coords with their index from the
  * latter sorted Coords.
  * @param q
  *   query Coord
  * @param ssrt
  *   subject SORTED Coords with index, we do not check
  *   if the ssrt is sorted or not in the function.
  * @param from
  *   [0, from) Coords from subject will be ignored. It
  *   is indepedent with index in subject.
  * @return
  *   Vector of subject coords and original index that
  *   overlap with the query Coord, with another index
  *   describing the position of subject. The latter
  *   index is mainly used for the next function
  *   [[findOvlp]] . The order of the subject coords
  *   are kept. The output may have zero element, use
  *   nonempty to check it.
  * TODO: use binary search instead
  */
def findItOvlp(q: Coord, ssrt: Vector[(Coord, Int)],
  from: Int): Vector[((Coord, Int), Int)] = {
  // zipWithIndex to record the position
  // on the current ssrt.
  ssrt.zipWithIndex
    .drop(from)
    .dropWhile(x => !isOvlp(x._1._1, q))
    .takeWhile(x => isOvlp(x._1._1, q))
}

/** Returns vector of objects depend on querie's
  * overlapping context
  *
  * Find Overlaped Coords in query and subject, then
  * perform a general function on the result. Same size
  * as qsrt will be returned if skipEmpty is false.
  *
  * @param qsrt
  *   Coord with index. The index may start from 0, or
  *   may be not ordered with arbitrary numbers. Here
  *   we keep the index since this function will often
  *   be used for The query Coords should be ordered by
  *   Coord (not index), which is checked inside the
  *   function. [[GenomeCoord]] objects, where we hope
  *   to perform overlap on each chrom and keep the
  *   input order finally. We will assign the index
  *   before performing the overlap function per
  *   chromosome.
  * @param ssrt
  *   Coord with index. The index involved here follows
  *   the same situation as discussed above. The
  *   subject Coords should be ordered by Coord (not
  *   index), which is checked, too.
  * @param f
  *   function used to get the result once the query
  *   has overlapping subjects The index here are the
  *   same as input.
  * @param g
  *   function used to get the result if no overlapping
  *   for a query. The second integer here is the
  *   scanned position in the given subject. Not the
  *   index from subject.
  * @param skipEmpty
  *   If no overlapping for a query, should we perform
  *   g or not. If false (not skip, default), we then
  *   have the same size of qsrt for output. If true
  *   (skip), only the queries with overlapping
  *   subjects have results.
  */
def findOvlp[S](qsrt: Vector[(Coord, Int)],
  ssrt: Vector[(Coord, Int)],
  f: ((Coord, Int), Vector[(Coord, Int)]) => S,
  g: ((Coord, Int), Int) => S,
  skipEmpty: Boolean = false): Vector[S] = {
  if (
    !isCoordSorted(qsrt.map(_._1)) &&
    !isCoordSorted(ssrt.map(_._1))
  ) {
    throw new RuntimeException("qsrt or ssrt is not sorted.")
  }
  val r = ListBuffer.empty[S]
  qsrt
    .foldLeft[Int](0)((si, q) => {
      val sovlp = findItOvlp(q._1, ssrt, from = si)
      if (sovlp.nonEmpty) {
        r.addOne(f(q, sovlp.map(_._1)))
        sovlp.head._2
      } else {
        if (!skipEmpty) {
          r.addOne(g(q, si))
        }
        si
      }
    }) // end of foldLeft
  r.toVector
}

/** Returns vector of boolean with the same size of
  * query.
  *
  * Check if query has overlapped subjects or not. Use
  * the [[findOvlp]] above with predefined f ang d,
  * skipEmpty is false.
  *
  * @param query
  *   sorted Coords
  * @param subject
  *   sorted Coords
  */
def isOverlapSorted(query: Coords,
  subject: Coords): Vector[Boolean] = {
  val f = (x: (Coord, Int), y: Vector[(Coord, Int)]) => true
  val g = (x: (Coord, Int), y: Int) => false
  findOvlp(query.zipWithIndex, subject.zipWithIndex, f, g,
      skipEmpty = false)
}

/** Returns query index and the corresponding
  * overlapped subject index.
  *
  * Only queries with overlaps returned, and the
  * subject index for each query is organized as
  * [FromSubjectIndex, ToSubjectIndex + 1).
  *
  * Finds Overlapped subjects for queries Use the
  * [[findOvlp]] above with predefined f and g,
  * skipEmpty is true.
  *
  * TODO: make reults as Vector[(Coord, Int)] to keep
  * consistent with other functions when using
  * zipWithIndex.
  *
  * @param query
  *   sorted Coords for mapping.
  * @param subject
  *   sorted Coords as pool for finding overlaps.
  */
def findOvlpOneChrSorted(query: Coords,
  subject: Coords): Vector[(Int, Coord)] = {
  val f = (x: (Coord, Int), y: Vector[(Coord, Int)]) =>
    (x._2, (y.head._2, y.last._2 + 1).asInstanceOf[Coord])
  val g = (x: (Coord, Int), y: Int) =>
    (x._2, (0, 0).asInstanceOf[Coord])
  findOvlp(query.zipWithIndex, subject.zipWithIndex, f, g,
      skipEmpty = true)
}

/** Returns the first overlapped subject for each
  * query. If no overlap found, (0, 0) Coords returned.
  *
  * The results is not unique and sorted, instead, it's
  * the same size as query.
  *
  * It calls [[findOvlp]] with predefined f and g with
  * skipEmpty as false.
  *
  * @param query
  *   sorted Coords for mapping
  * @param subject
  *   sorted Coords as pool for finding overlaps.
  */
def findSubject(query: Coords, subject: Coords): Coords = {
  val f = (x: (Coord, Int), y: Vector[(Coord, Int)]) => y.head._1
  val g = (x: (Coord, Int), y: Int) => (0, 0).asInstanceOf[Coord]
  findOvlp(query.zipWithIndex, subject.zipWithIndex, f, g,
      skipEmpty = false)
}

/** If gene's strand is "+", then
  *       max(0, startFrom - upStream)
  *        |      min(startFrom + downStream, chromSize)
  *        | promoter >> |
  * 5' ---------------------> 3'
  *            |
  *        startFrom is TSS
  *
  * If gene's strand is "-", then
  *       max(0, endTo - downStream)
  *        |      min(endTo + upStream, chromSize)
  *        | << promoter |
  * 5' ---------------------> 3'
  *               |
  *        endTo is TSS
  * If strand is others, then will be treated as "+".
  */
trait CoordOpt {
  def getUpStream(x: Coord, upStream: Int): Int
  def getDownStream(x: Coord, downStream: Int): Int
}

object GeneBody3to5 extends CoordOpt {
  override def getDownStream(x: Coord, downStream: Int): Int = {
    (x.startFrom - downStream).max(0)
  }
  override def getUpStream(x: Coord, upStream: Int): Int = {
    x.endTo + upStream
  }
}

object GeneBody5to3 extends CoordOpt {
  override def getDownStream(x: Coord, downStream: Int): Int = {
    x.endTo + downStream
  }
  override def getUpStream(x: Coord, upStream: Int): Int = {
    (x.startFrom - upStream).max(0)
  }
}

object Promoter3to5 extends CoordOpt {
  override def getDownStream(x: Coord, downStream: Int): Int = {
    (x.endTo - downStream).max(0)
  }
  override def getUpStream(x: Coord, upStream: Int): Int = {
    x.endTo + upStream
  }
}

object Promoter5to3 extends CoordOpt {
  override def getDownStream(x: Coord, downStream: Int): Int = {
    x.startFrom + downStream
  }
  override def getUpStream(x: Coord, upStream: Int): Int = {
    (x.startFrom - upStream).max(0)
  }
}

/** Returns the closed coords (one or more) for a given
  * coord.
  *
  *   1. Empty Vector if all coords are beyond the
  *      range.
  *   2. Closeness is defined by coords distance
  *      [[d()]].
  *   3. All overlapped coords (if had) will be
  *      returned, else only closest one.
  * @param x
  *   the coord for closed coords
  * @param y
  *   the pool of coords for selection
  * @param within
  *   upstream and downstream of x's two sides
  */
def getClosedCoord(x: Coord, y: Coords,
  within: Int = 1e6.toInt): Coords = {
  val leftMost: Int  = (x.startFrom - within).max(0)
  val rightMost: Int = x.endTo + within
  val p              = y
    .dropWhile(t => {
      (t.endTo < leftMost) || (t.startFrom > rightMost)
    })
    .sorted(using coordOrd)
  if (p.isEmpty) {
    Vector[Coord]()
  } else {
    val r     = p.map(t => (t, d(x, t))).sortBy(_._2)
    val index = r.indexWhere(x => x._2 > 0)
    if (index < 0) {
      r.map(x => x._1)
    } else {
      r.take(index.max(1)).map(x => x._1)
    }
  }
}
