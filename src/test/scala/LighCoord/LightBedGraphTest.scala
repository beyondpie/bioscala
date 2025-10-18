import org.scalatest._
import org.scalatest.funsuite.AnyFunSuite
import os._
import scala.io.Source
import bioscala.LightCoord.BedGraph._
import SZUtils.{str2path, path2str}
import scala.util.Try
import scala.util.Success
import scala.util.Failure
import bioscala.LightCoord.GenomeCoord.GenomeCoord
import java.io.FileNotFoundException
import org.broad.igv.bbfile.BBFileReader
import scala.collection.mutable.ListBuffer

class LightBedGraphTest extends AnyFunSuite {
  val bedgraphfnm = getClass.getResource("/test.bedgraph").getPath
  val chromsizefnm = getClass.getResource("/test.chromSize.txt").getPath
  val bwfnm = getClass.getResource("/test.bw").getPath

  val g1: GenomeCoord = (
    chrom = "chr1", coord = (startFrom = 100, endTo = 101),
    strand = ".")
  val g2: GenomeCoord = (
    chrom = "chr1", coord = (startFrom = 200, endTo = 400),
    strand = "."
  )
  val g3: GenomeCoord = (
    chrom = "chr3", coord = (startFrom = 0, endTo = 100),
    strand = "."
  )

  test("Test getBigWigReader") {
    val r = getBigWigReader(bwfnm)
    val y1 = getWigValue(r, g1, 0.0, false)
    println(y1)
    assert(y1.head.s > 0.0)

    val y3 = getWigValue(r, g3, -1.0, false)
    println(y3)
    assert(y3.head.s < 0.0)
  }

  test("Test getbwOneRegion") {
    val r = getBigWigReader(bwfnm)
    val y2 = getbwOneRegion(r, g2, 0.0, true, false)
    println(s"weight by region size: $y2")
    assert(y2 == 0.45)
    val y22 = getbwOneRegion(r, g2, 0.0, false, false)
    println(s"weight by each ovlp region size: $y22")
    assert(y22 == 0.75)
  }

}

