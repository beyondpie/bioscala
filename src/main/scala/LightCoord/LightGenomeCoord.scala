package bioscala.LightCoord.GenomeCoord

import Genome.GenomeOpts.chromOrd
import bioscala.LightCoord._
import bioscala.LightCoord.coordOrd
import os._


type GenomeCoord  = (chrom: String, coord: Coord, strand: String)
type GenomeCoords = Vector[GenomeCoord]

given genomeCoordIgnoreStrandOrd: Ordering[GenomeCoord] with {
  def compare(x: GenomeCoord, y: GenomeCoord): Int = {
    if (x.chrom == y.chrom) {
      coordOrd.compare(x.coord, y.coord)
    } else {
      chromOrd.compare(x.chrom, y.chrom)
    }
  }
}

// TODO: use extension methods
def toStringGenomeCoord(g: GenomeCoord,
  useStrand: Boolean = false): String = {
  val base = s"${g.chrom}:${toStringCoord(g.coord)}"
  if (useStrand) {
    s"$base:${g.strand}"
  } else {
    base
  }
}

// TODO: use extension methods
def mkStringGenomeCoord(x: GenomeCoord, sep: String = "\t",
  useStrand: Boolean = true): String = {
  val base =
    s"${x.chrom}$sep${x.coord.startFrom}$sep${x.coord.endTo}"
  if (useStrand) {
    s"$base$sep${x.strand}"
  } else {
    base
  }
}


def fromBed(bed: String): GenomeCoords = {
  os.read.lines.stream(os.Path(bed))
    .map(x => x.strip().split("\t"))
    .filter(x => x.nonEmpty)
    .filter(x => x.length >= 3)
    .map(x => (chrom = x.head,
      coord = (x(1).toInt, x(2).toInt),
      strand = "."
    )).toVector
}

/** Returns processed overlaping results for
  * GenomeCoords, which is the same size and order with
  * the query.
  *
  * A general function stands on
  * [[bioscala.LightCoord.findOvlp]] with skipEmpty as
  * false. for each chromosome in query.
  *
  * @param query
  *   Genome Coords without sorting, from which we map
  *   to a given pool of Genome Coords to get the
  *   overlapped results.
  * @param subject
  *   Genome Coords without sorting, to which we map
  *   our candidates to.
  * @param f0
  *   see details at [[bioscala.LightCoord.findOvlp]]
  * @param f1
  *   used to process the results from f0.
  * @param g0
  *   see details at [[bioscala.LightCoord.findOvlp]]
  * @param g1
  *   used to process the reulsts if no chromosome
  *   matched in the subject
  */
def findOvlpGenomeCoordIgnoreStrand[S, T](
  query: GenomeCoords,
  subject: GenomeCoords,
  f0: ((Coord, Int), Vector[(Coord, Int)]) => S,
  f1: ((GenomeCoord, Int), S) => T,
  g0: ((Coord, Int), Int) => S,
  g1: (GenomeCoord, Int) => T
): Vector[T] = {
  val ss = subject.zipWithIndex.groupBy(_._1.chrom)
  query.zipWithIndex
    .groupBy(_._1.chrom)
    .map((chrom, q) => {
      if (ss.keySet.contains(chrom)) {
        val qq = q.sortBy(_._1.coord)(using coordOrd)
        val sord = ss(chrom)
          .map(x => (x._1.coord, x._2))
          .sortBy(x => x._1)(using coordOrd)
        val x = findOvlp(qq.map(x => (x._1.coord, x._2)), sord,
            f0, g0, skipEmpty = false)
        qq.zip(x).map((g, y) => (g._2, f1(g, y)))
      } else {
        q.map(x => (x._2, g1(x._1, x._2)))
      }
    })
    .flatten
    .toVector
    .sortBy(_._1)
    .map(_._2)
}

/** Returns vector of boolean to check if there is any
  * GenomeCoord overlapped with a given element in
  * query. Keep both the size and order for query.
  *
  * Uses[[findOvlpGenomeCoordIgnoreStrand]] directly
  * with predefiend functions.
  *
  * @param query
  *   Genome Coords to check without ordering
  * @param subject
  *   Genome Coords to be mapped to without ordering
  */
def isOverlapIgnoreStrand(
  query: GenomeCoords,
  subject: GenomeCoords
): Vector[Boolean] = {
  val f0 = (x: (Coord, Int), y: Vector[(Coord, Int)]) => true
  val g0 = (x: (Coord, Int), y: Int) => false
  val f1 = (x: (GenomeCoord, Int), y: Boolean) => y
  val g1 = (x: GenomeCoord, y: Int) => false
  findOvlpGenomeCoordIgnoreStrand[Boolean, Boolean](
      query,
      subject,
      f0,
      f1,
      g0,
      g1
  )
}

/** Returns vector of GenomeCoord to get the first
  * overlapped Genome Coord from subject for a given
  * element in query. Keep both the size and order for
  * query. If no overlapped Genome Coords from the
  * subject, a (0, 0) Coord returned on the same
  * chromosome.
  *
  * Use [[findOvlpGenomeCoordIgnoreStrand]] directly
  * with predefined functions.
  *
  * @param query
  *   Geneome Coords from which we map without ordering
  * @param subject
  *   Genom Coords to be mapped to without ordering
  */
def findSubjectOvlpGenomeCoordIgnoreStrand(
  query: GenomeCoords,
  subject: GenomeCoords
): GenomeCoords = {
  val f0 = (x: (Coord, Int), y: Vector[(Coord, Int)]) => y.head._1
  val f1 = (x: (GenomeCoord, Int), y: Coord) => {
    x._1.toTuple.copy(_2 = y).asInstanceOf[GenomeCoord]
  }
  val g0 = (x: (Coord, Int), y: Int) => (startFrom = 0, endTo = 0)
  val g1 = (x: GenomeCoord, y: Int) => {
    x.toTuple.copy(_2 = (startFrom = 0, endTo = 0))
  }
  findOvlpGenomeCoordIgnoreStrand(
      query,
      subject,
      f0,
      f1,
      g0,
      g1
  )
}

