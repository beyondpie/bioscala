package bioscala.LightCoord

import scala.collection.mutable.ListBuffer
import SZUtils.ifelse

/**
 * Genome element range in a given chrom.
 * - Only start and end position, no chrom info.
 * - Always from 5' to 3', i.e., startFrom <= endTo
 * - Follow Bed format, 0-index and range will be left closed and right open.
 * - both start and end positions should be no smaller than zero.
 */
type Coord  = (startFrom: Int, endTo: Int)
type Coords = Vector[Coord]

/**
 * Get the center position of a coord.
 * Here center is startFrom + |__ (endTo - startFrom) / 2 __|
 * The latter part will be integer part of the division.
 * @param x single coord
 * @return
 */
def getCenter(x: Coord): Coord = {
  val c: Int = (x.endTo - x.startFrom) / 2
  (startFrom = c, endTo = c + 1)
}

/**
 * If gene's strand is "+", then
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

def toStringCoord(g: Coord): String = {
  s"${g.startFrom.toInt}-${g.endTo.toInt}"
}

def mkStringCoord(x: Coord, sep: String = "\t"): String = {
  s"${x.startFrom}${sep}${x.endTo}"
}

// TODO: chrom2int should be moved to chrom or genome
def chrom2int(x: String): Int = {
  val xchr = x.toLowerCase().replace("chr", "")
  if (xchr == "X") {
    100
  } else if (xchr == "Y") {
    101
  } else {
    xchr.toInt
  }
}

given chromOrd: Ordering[String] with {
  def compare(x: String, y: String): Int = {
    chrom2int(x).compare(chrom2int(y))
  }
}

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

/**
 * Check if a vector of genomic coordinates are well sorted or not Ref:
 * StackOverflow about check-whether-a-collection-is-ordered
 * @param x
 * @return
 */
def isCoordSorted(x: Iterable[Coord]): Boolean = {
  x.view.zip(x.tail).forall(((x, y) => coordOrd.lteq(x, y)))
}

/**
 * Intersect between two coords.
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
      Some((startFrom = y.startFrom, endTo = x.endTo.min(y.endTo)))
    }
  } else {
    if (y.endTo <= x.startFrom) {
      None
    } else {
      Some(startFrom = x.startFrom, endTo = y.endTo.min(x.endTo))
    }
  }
}

/**
 * Find Overlap between two ranges in the same chromosomes. The input
 * query and subject must be Sorted.
 *   - strand is ignored.
 *   - as long as they share at least one base pair (bed format)
 *
 * Implementation details: We find the first element position in subject
 * that overlaps with a given query by using `dropWhile`; then we keep
 * the elements untill the elements have no overlaps with the query term
 * by using `takeWhile`. If elements left, we then keep the coordinates
 * of the subject one.
 *
 * @param query
 *   the coords we want to check if they have overlaps on another pool
 *   of coord
 * @param subject
 *   the pool of coords for query coords to map to
 * @return
 *   QueryIndex and a range of Subject Index overlapping with the query.
 *   1. Only the index of query overlapped with subject will be
 *      returned.
 *   2. QueryIndex and subjectIndex, start from 0. 3. The
 *      overlapped genomic coordinates from Subject are recorded as
 *      [startFrom, endTo + 1), i.e., left closed and right open 4. Here
 *      we use Coord as the range of subject Index, which is not the
 *      real coord. TODO: use (Int, Int) instead of Coord to explicitly
 *      express what we are doing.
 */
def findOvlpOneChrSorted(query: Coords,
  subject: Coords): Vector[(Int, Coord)] = {
  if (!isCoordSorted(query)) {
    throw new RuntimeException("query is not sorted.")
  }
  if (!isCoordSorted(subject)) {
    throw new RuntimeException("subject is not sorted.")
  }
  val r = ListBuffer.empty[(Int, Coord)]
  query.zipWithIndex.foldLeft[Int](0)((si, q) => {
    val sovlp =
      subject.zipWithIndex
        .dropWhile(((x, id) => {
            (id < si) ||
            (x.endTo < q._1.startFrom) || (x.startFrom >= q._1.endTo)
          }))
        .takeWhile(((x, id) => {
            (x.endTo >= q._1.startFrom) && (x.startFrom < q._1.endTo)
          }))
    if (sovlp.nonEmpty) {
      r.addOne((q._2, (sovlp.head._2, sovlp.last._2 + 1)))
      sovlp.head._2
    } else {
      si
    }
  })
  r.toVector
}

// TODO: this function has lots of copied codes from the function above,
// which we can provide one unified API.
def isOverlapSorted(query: Coords, subject: Coords): Vector[Boolean] = {
  if (!isCoordSorted(query)) {
    throw new RuntimeException("query is not sorted.")
  }
  if (!isCoordSorted(subject)) {
    throw new RuntimeException("subject is not sorted.")
  }
  val r = ListBuffer.empty[Boolean]
  query.foldLeft[Int](0)((si, q) => {
    val sovlp = subject.zipWithIndex
      .dropWhile(((x, id) => {
          (id < si) || (x.endTo < q.startFrom) || (x.startFrom) >= q.endTo
        }))
      .takeWhile(((x, id) => {
          (x.endTo >= q.startFrom) && (x.startFrom < q.endTo)
        }))
    if (sovlp.nonEmpty) {
      r.addOne(true)
      sovlp.head._2
    } else {
      r.addOne(false)
      si
    }
  })
  r.toVector
}

// General ovlp method to test

def isOvlp(x: Coord, y: Coord): Boolean = {
  ifelse(x.endTo < y.startFrom || x.startFrom >= y.endTo, false, true)
}

def findOvlp(q: Coord, ssrt: Vector[(Coord, Int)],
  from: Int): Vector[(Coord, Int)] = {
  ssrt
    .drop(from)
    .dropWhile(x => !isOvlp(x._1, q))
    .takeWhile(x => isOvlp(x._1, q))
}

def findOvlp[S](qsrt: Vector[(Coord, Int)], ssrt: Vector[(Coord, Int)],
  f: ((Coord, Int), Vector[(Coord, Int)]) => S,
  g: ((Coord, Int), Int) => S,
  skipEmpty: Boolean = false): Vector[S] = {
  if (!isCoordSorted(qsrt.map(_._1)) &&
    !isCoordSorted(ssrt.map(_._1))) {
    throw new RuntimeException("qsrt or ssrt is not sorted.")
  }
  val r      = ListBuffer.empty[S]
  qsrt
    .foldLeft[Int](0)((si, q) => {
      val sovlp: Vector[(Coord, Int)] =
        findOvlp(q._1, ssrt, from = si)
      if (sovlp.nonEmpty) {
        r.addOne(f(q, sovlp))
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

def isOverlapSorted2(query: Coords,
  subject: Coords): Vector[Boolean] = {
  val f = (x: (Coord, Int), y: Vector[(Coord, Int)]) => true
  val g = (x: (Coord, Int), y: Int) => false
  findOvlp(query.zipWithIndex, subject.zipWithIndex, f, g, skipEmpty = false)
}

def findOvlpOneChrSorted2(query: Coords,
  subject: Coords): Vector[(Int, Coord)] = {
  val f = (x: (Coord, Int), y: Vector[(Coord, Int)]) =>
    (x._2, (y.head._2, y.last._2 + 1).asInstanceOf[Coord])
  val g = (x: (Coord, Int), y: Int) =>
    (x._2, (0, 0).asInstanceOf[Coord])
  findOvlp(query.zipWithIndex, subject.zipWithIndex, f, g, skipEmpty = true)
}

def findSubject(query: Coords, subject: Coords): Coords = {
  val f = (x: (Coord, Int), y: Vector[(Coord, Int)]) => y.head._1
  val g = (x: (Coord, Int), y: Int) => (0, 0).asInstanceOf[Coord]
  findOvlp(query.zipWithIndex, subject.zipWithIndex, f, g, skipEmpty = true)
}
