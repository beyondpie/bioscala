package GenomicCoordinate
// will be stable on Scala version >= 3.7.0
import scala.language.experimental.namedTuples
import scala.collection.mutable.ListBuffer

type Coord = (startFrom: Int, endTo: Int)


/**
  * Check if a vector of genomic coordinates are well sorted or not
  * Ref: StackOverflow about check-whether-a-collection-is-ordered
  * @param x
  * @return
  */
def isGenomicCoordSorted(x: Iterable[Coord]):Boolean = {
  x.view.zip(x.tail).toStream.forall( (x, y) => {
    x.startFrom <= y.startFrom 
  })
}


/**
  * Find Overlap between two ranges in the same chromosomes.
  * The input query and subject must be Sorted.
  * - strand is ignored.
  * - as long as they share at least one base pair (bed format)
  *
  * Implementation details:
  * We find the first element position in subject that overlaps with
  * a given query by using `dropWhile`; then we keep the elements
  * untill the elements have no overlaps with the query term by using
  * `takeWhile`. If elements left, we then keep the coordinates of the
  * subject ones.
  * 
  * @param a
  * @param b
  * @return
  * 1. Only queryIndex with overlapped subjectIndex will be returned.
  * 2. QueryIndex and subjectIndex, start from 0.
  * 3. The overlapped genomic coordinates from Subject are recorded as
  *    [startId, endId + 1), i.e., left closed and right open
  * 
  */
def findOvlpOneChrSorted(query: Vector[Coord],
  subject: Vector[Coord]): Vector[(Int, Coord)] = {
  if (!isGenomicCoordSorted(query)) {
    throw new RuntimeException("query is not sorted.")
  }
  if (!isGenomicCoordSorted(subject)) {
    throw new RuntimeException("subject is not sorted.")
  }
  val r = ListBuffer.empty[(Int, Coord)]
  query.zipWithIndex.foldLeft[Int](0)((si, q) => {
    val sovlp =
      subject.zipWithIndex
        .dropWhile((x, id) => {
          (id < si) ||
            (x.endTo < q._1.startFrom) || (x.startFrom >= q._1.endTo)
        })
        .takeWhile((x, id) => {
          (x.endTo >= q._1.startFrom) && (x.startFrom < q._1.endTo)
        })
    if (sovlp.length >= 1) {
      r.addOne((q._2, (sovlp.head._2, sovlp.last._2 + 1)))
      sovlp.head._2
    } else {
      si
    }
  })
  r.toVector
}


/**
  * Intersect between two coords.
  *
  * @param x
  * @param y
  * @return
  */
def intersect(x: Coord, y: Coord): Option[Coord] = {
  if (x.startFrom <= y.startFrom ) {
    if (x.endTo <= y.startFrom) {
      None
    }
    else {
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
