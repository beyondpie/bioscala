package bioscala.LightCoord

type GenomeCoord  = (chrom: String, coord: Coord, strand: String)
type GenomeCoords = Vector[GenomeCoord]

def toStringGenomeCoord(g: GenomeCoord,
  useStrand: Boolean = false): String = {
  val base = s"${g.chrom}:${toStringCoord(g.coord)}"
  if (useStrand) {
    s"$base:${g.strand}"
  } else {
    base
  }
}

def mkStringGenomeCoord(x: GenomeCoord, sep: String = "\t",
  useStrand: Boolean = true): String = {
  val base = s"${x.chrom}$sep${x.coord.startFrom}$sep${x.coord.endTo}"
  if (useStrand) {
    s"$base$sep${x.strand}"
  } else {
    base
  }
}

given genomeCoordIgnoreStrandOrd: Ordering[GenomeCoord] with {
  def compare(x: GenomeCoord, y: GenomeCoord): Int = {
    if (x.chrom == y.chrom) {
      coordOrd.compare(x.coord, y.coord)
    } else {
      chromOrd.compare(x.chrom, y.chrom)
    }
  }
}

def isOverlapIgnoreStrand(query: GenomeCoords,
  subject: GenomeCoords): Vector[Boolean] = {
  val ss: Map[String, GenomeCoords] = subject.groupBy(_.chrom)
  query.zipWithIndex
    .groupBy(_._1.chrom)
    .flatMap(
        ((chrom, q) => {
            if (ss.keySet.contains(chrom)) {
              val qq   = q.sortBy(_._1.coord)(using coordOrd)
              val sord = ss(chrom).map(_.coord).sorted(using coordOrd)
              qq.map(_._2)
                .zip(isOverlapSorted(qq.map(_._1.coord), sord))
            } else {
              q.map(x => (x._2, false))
            }
          }))
    .toVector
    .sortBy(_._1)
    .map(_._2)
}

// a general implementation of findOvlp
def findOvlpGenomeCoordIgnoreStrand[S, T](
  query: GenomeCoords,
  subject: GenomeCoords,
  f0: ((Coord, Int), Vector[(Coord, Int)]) => S,
  f1: ((GenomeCoord, Int), S) => T,
  g0: ((Coord, Int), Int) => S,
  g1: (GenomeCoord, Int) => T
): Vector[T] = {
  val ss: Map[String, Vector[(GenomeCoord, Int)]] =
    subject.zipWithIndex.groupBy(_._1.chrom)
  query.zipWithIndex
    .groupBy(_._1.chrom)
    .flatMap(
        ((chrom, q) => {
            if (ss.keySet.contains(chrom)) {
              val qq = q
                .sortBy(_._1.coord)(using coordOrd)
              val sord = ss(chrom)
                .map(x => (x._1.coord, x._2))
                .sortBy(x => x._1)(using coordOrd)
              val x0: Vector[S] =
                findOvlp(qq.map(x => (x._1.coord, x._2)), sord, f0, g0,
                    skipEmpty = false)
              qq.zip(x0).map((q, x) => (f1(q, x), q._2))
            } else {
              q.map(x => g1(x._1, x._2)).zip(q.map(_._2))
            }
          }))
    .toVector
    .sortBy(_._2)
    .map(_._1)
}

def isOverlapGenomeCoordIgnoreStrand(
  query: GenomeCoords,
  subject: GenomeCoords
): Vector[Boolean] = {
  val f0 = (x: (Coord, Int), y: Vector[(Coord, Int)]) => true
  val f1 = (x: (GenomeCoord, Int), y: Boolean) => true
  val g0 = (x: (Coord, Int), y: Int) => false
  val g1 = (x: GenomeCoord, y: Int) => false
  findOvlpGenomeCoordIgnoreStrand(
      query,
      subject,
      f0,
      f1,
      g0,
      g1
  )
}

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
