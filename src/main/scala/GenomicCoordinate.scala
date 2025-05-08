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
def findOvlpOneChrSorted(query: Vector[(Int, Int)], subject: Vector[(Int, Int)]): Vector[(Int, (Int, Int))] = {
  if (!isGenomicCoordSorted(query)) {
    throw new RuntimeException("query is not sorted.")
  }
  if (!isGenomicCoordSorted(subject)) {
    throw new RuntimeException("subject is not sorted.")
  }
  val r = ListBuffer.empty[(Int, (Int, Int))]
  query.zipWithIndex.foldLeft[Int](0)((si, q) => {
    val sovlp =
      subject.zipWithIndex
        .dropWhile((x, id) => {
          (id < si) || (x._2 < q._1._1) || (x._1 >= q._1._2)
        })
        .takeWhile((x, id) => {
          (x._2 >= q._1._1) && (x._1 < q._1._2)
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
