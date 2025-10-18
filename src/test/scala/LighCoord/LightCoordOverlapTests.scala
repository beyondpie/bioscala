import org.scalatest._
import org.scalatest.funsuite.AnyFunSuite
import bioscala.LightCoord.isOverlapSorted as isOverlapSorted2
import bioscala.LightCoord.{Coord, Coords}
import bioscala.LightCoord.coordOrd
import bioscala.LightCoord.isCoordSorted
import SZUtils.ifelse
import bioscala.LightCoord.findOvlpOneChrSorted as findOvlpOneChrSorted2
import scala.io.Source
import bioscala.LightCoord.findSubject
import scala.collection.mutable.ListBuffer

class CoordFindOvlpSuite extends AnyFunSuite {
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

  def isOverlapSorted(query: Coords,
    subject: Coords): Vector[Boolean] = {
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

  def getCoordFromResource(fnm: String,
    chr: String = "chr1"): Coords = {
    Source
      .fromResource(fnm)
      .getLines
      .drop(1)
      .map(x => x.strip().split("\t"))
      .filter(x => x(0) == chr)
      .map(x => (startFrom = x(1).toInt, endTo = x(2).toInt))
      .toVector
      .sorted(using coordOrd)
  }

  for (chr <- List("chr1", "chr2", "chr10")) {
    val c1 = getCoordFromResource(
        "001_CLA_EPd_CTX_Car3_Glut-H3K4me1-Male.reproPeak", chr)
    val c2 = getCoordFromResource(
        "002_IT_EP_CLA_Glut-H3K4me1-Male.reproPeak",
        chr
    )
    test(s"Coord isOverlap for $chr") {
      val a1 = isOverlapSorted(c1, c2).map(x => ifelse(x, 1, 0))
      val a2 = isOverlapSorted2(c1, c2).map(x => ifelse(x, 1, 0))
      println(s"In $chr, ${a1.sum} a1 and ${a2.sum} a2")
      assert(a1.zip(a2).map((x, y) => x - y).sum == 0)
    }

    test(s"Coord findOverlap for $chr") {
      val a1 = findOvlpOneChrSorted(c1, c2)
      val a2 = findOvlpOneChrSorted2(c1, c2)
      // println(s"Im $chr, ${a1.length} a1 and ${a2.length} a2.")
      assert(a1.zip(a2).map((x, y) => (x._1 - y._1)).sum == 0)
      assert(a1
            .zip(a2)
            .forall((x, y) =>
              x._2.startFrom == y._2.startFrom &&
                x._2.endTo == y._2.endTo))
    } // ned of test
    test(s"Coord findSubject for $chr") {
      val a1 = findSubject(c1, c2)
        .distinct.sorted(using coordOrd)
      val a2 = findOvlpOneChrSorted2(c1, c2)
      println(c1(a2.head._1))
      println(c2(a2.head._2.startFrom))
      println(c2(a2.head._2.endTo - 1))
      println(a1.length)
      println(a1.head)
      assert(a1.length > 0)
    }

  } // end of for
}
