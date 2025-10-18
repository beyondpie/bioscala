import org.scalatest.funsuite.AnyFunSuite
import bioscala.LightCoord.coordOrd
import bioscala.LightCoord.isOverlapSorted
import bioscala.LightCoord.GenomeCoord._
import bioscala.LightCoord.GenomeCoord.isOverlapIgnoreStrand as isOverlapIgnoreStrand2
import SZUtils.ifelse
import scala.io.Source

class GenomeCoordFindOvlpSuite extends AnyFunSuite {
  def isOverlapIgnoreStrand(query: GenomeCoords,
    subject: GenomeCoords): Vector[Boolean] = {
    val ss: Map[String, GenomeCoords] = subject.groupBy(_.chrom)
    query.zipWithIndex
      .groupBy(_._1.chrom)
      .flatMap((chrom, q) => {
        if (ss.keySet.contains(chrom)) {
          val qq   = q.sortBy(_._1.coord)(using coordOrd)
          val sord = ss(chrom).map(_.coord).sorted(using coordOrd)
          qq.map(_._2)
            .zip(isOverlapSorted(qq.map(_._1.coord), sord))
        } else {
          q.map(x => (x._2, false))
        }
      })
      .toVector
      .sortBy(_._1)
      .map(_._2)
  }

  def getGenomeCoordFromResource(fnm: String,
    chrs: Vector[String]): GenomeCoords = {
    Source
      .fromResource(fnm)
      .getLines
      .drop(1)
      .map(x => x.strip().split("\t"))
      .filter(x => chrs.contains(x(0)))
      .map(x =>
        (chrom = x(0),
            coord = (startFrom = x(1).toInt, endTo = x(2).toInt),
            strand = "."))
      .toVector
  }

  // val chrs: Vector[String] = Vector("chr1", "chr2", "chr10")
  val chrs: Vector[Vector[String]] = Vector(
      Vector("chr1"),
      Vector("chr2"),
      Vector("chr10"),
      Vector("chr1", "chr2"),
      Vector("chr1", "chr10"),
      Vector("chr1", "chr2", "chr10"),
      Vector("chr2", "chr10"),
      Vector("chr10", "chr2")
  )
  for (chr <- chrs) {
    val c1: GenomeCoords = getGenomeCoordFromResource(
        fnm = "001_CLA_EPd_CTX_Car3_Glut-H3K4me1-Male.reproPeak",
        chr)
    val c2: GenomeCoords = getGenomeCoordFromResource(
        fnm = "002_IT_EP_CLA_Glut-H3K4me1-Male.reproPeak", chr)
    val chrstr = chr.mkString(",")
    test(s"GenomeCoord isOverlapIgnoreStrand for ${chrstr}") {
      val a1: Vector[Int] =
        isOverlapIgnoreStrand(c1, c2).map(x => ifelse(x, 1, 0))
      val a2: Vector[Int] =
        isOverlapIgnoreStrand2(c1, c2).map(x => ifelse(x, 1, 0))
      println(s"a1 size ${a1.length} and a2 size ${a2.length}")
      println(s"${a1.sum} a1 and ${a2.sum} a2")
      assert(a1.zip(a2).map((x, y) => (x - y).abs).sum == 0)
    }
  }

  for (chr <- chrs) {
    val c1: GenomeCoords = getGenomeCoordFromResource(
        fnm = "001_CLA_EPd_CTX_Car3_Glut-H3K4me1-Male.reproPeak",
        chr)
    val c2: GenomeCoords = getGenomeCoordFromResource(
        fnm = "002_IT_EP_CLA_Glut-H3K4me1-Male.reproPeak", chr)
    val chrstr = chr.mkString(",")
    test(s"GenomeCoord findSubjectOvlp for ${chrstr}") {
      val a2: Vector[Int] =
        isOverlapIgnoreStrand2(c1, c2).map(x => ifelse(x, 1, 0))
      val r  = findSubjectOvlpGenomeCoordIgnoreStrand(c1, c2)
      val p = c1.zip(r).filter((x, y) => y.coord.endTo != 0)
        .map((x, y) => toStringGenomeCoord(x) + "ovlp" + toStringGenomeCoord(y))
        .take(5)
      println(p)

      val np = c1.zip(r).filter((x, y) => y.coord.endTo == 0)
        .map((x, y) => toStringGenomeCoord(x) + "noovlp" + toStringGenomeCoord(y))
        .take(5)
      println(np)

      val a3 = r.map(x => ifelse(x.coord.endTo != 0, 1, 0))
      assert(a2.zip(a3).map((x, y) => (x - y).abs).sum == 0)
    }
  }
}
