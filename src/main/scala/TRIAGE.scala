package bioscala.TRIAGE

import Genome.MouseGenome.chr2size
import SZUtils.ifelse
import bioscala.LightCoord.Bed.LightBed
import bioscala.LightCoord.Bed.LightBedElement
import bioscala.LightCoord.Coord
import bioscala.LightCoord.Coords
import bioscala.LightCoord.GenCode.GFF3
import bioscala.LightCoord.GenCode.aroundGene
import bioscala.LightCoord.GeneBody3to5
import bioscala.LightCoord.GeneBody5to3
import bioscala.LightCoord.GenomeCoord.GenomeCoord
import bioscala.LightCoord.GenomeCoord.findOvlpGenomeCoordIgnoreStrand
import bioscala.LightCoord.GenomeCoord.genomeCoordIgnoreStrandOrd
import bioscala.LightCoord.coordOrd
import bioscala.LightCoord.findOvlpOneChrSorted
import bioscala.LightCoord.getCenter
import smile.math.MathEx.mean
import smile.math.MathEx.stdev

import math.*
import bioscala.LightCoord.GenCode.GFF3Element

/** Group Bed elements into chromosome majored, and
  * sort the coords within each chromosome.
  *
  * @param d
  *   vector of LightBedElements
  * @return
  */
def groupSortLightBed(d: LightBed): Map[String, LightBed] = {
  d.groupBy(_.g.chrom)
    .map((chr, v) => (chr, v.sortBy(_.g.coord)(using coordOrd)))
}

/** Assign Broadest domain to gene bodies. Ref: TRIAGE,
  * https://linkinghub.elsevier.com/retrieve/pii/S2405471220304191
  *
  * For each gene body, if broad domain found, will
  * assign the broadest domain as genome coord,
  * namefn(gene) as name, raw broad domain score
  * (domain width) as score. If no broad domain found,
  * will assign extended gene body as genome coord,
  * defaultScore as score.
  *
  *   1. Extend the gene body based on upStream to TSS,
  *      and downStream to TES as proximal region, see
  *      [[aroundGene]].
  *   2. Map the center / any position of broad domain
  *      to the proximal region.
  *   3. Only broadest domain used for each gene.
  *
  * In TRIGAGE, it seems that they firslty use domain
  * center for each gene's overlapping, Then later they
  * consider all the domains having at least one gene
  * mapped above, and check if the domains also overlap
  * with other genes though no center overlapped.
  *
  * @param x
  *   genes
  * @param d
  *   broad domain regions for mapping or overlapping
  *   with genes they can be centers of ranges (TRIGAGE
  *   mainly uses this) or anything else.
  * @param rd
  *   raw broad domain regions, used for assigning the
  *   genome coord once overlapping happens.
  * @param namefn
  *   to assign name for once the overlapping happens.
  * @param defaultScore
  *   to assign score if no overlapping found.
  */
def assignBroadDomain2Gene(genes: GFF3, d: LightBed, rd: LightBed,
  namefn: GFF3Element => String,
  defaultScore: Double = -1.0): LightBed = {
  val query   = genes.map(x => x.g)
  val subject = d.map(d => d.g)
  val f0 = (x: (Coord, Int), y: Vector[(Coord, Int)]) => {
    val yid = y.maxBy((coord, index) => d(index).score)._2
    (
        g = rd(yid).g, // use raw broad domain region
        name = namefn(genes(x._2)),
        score = d(yid).score
    )
  }
  val f1 = (x: (GenomeCoord, Int), y: LightBedElement) => y
  val g0 = (x: (Coord, Int), y: Int) => {
    (
        g = genes(x._2).g,
        name = namefn(genes(x._2)),
        score = defaultScore
    )
  }
  val g1 = (x: GenomeCoord, xid: Int) => {
    (
        g = genes(xid).g,
        name = namefn(genes(xid)),
        score = defaultScore
    )
  }
  findOvlpGenomeCoordIgnoreStrand(query, subject, f0, f1, g0, g1)
}

/** Get all the unique LightBedElements. The elements
  * will be ordered by genomeCoordIgnoreStrandOrd
  * @param x
  *   matrix of LightBedElement
  * @param defaultScore
  *   filter element by it (not included) before
  *   listing all of them
  * @return
  */
def listAllUniqueElement(x: Vector[LightBed],
  defaultScore: Double = 0.0): LightBed = {
  x.flatMap(y => y.filter(g => g.score > defaultScore))
    .distinct
    .sortBy(_.g)(using genomeCoordIgnoreStrandOrd)
}

/** Get the top domain width from the mapped domains
  *
  * Ref:
  * https://linkinghub.elsevier.com/retrieve/pii/S2405471220304191
  *
  * @param x
  *   assigned domains organized as subclass by genes
  * @param qtop
  *   top ratio, which should be determined by the
  *   coverages of variably expressed TFs (VETFs) or
  *   other genes using the top domains. Follow the
  *   reference, we use 0.05 as default.
  * @param defaultScore
  *   used for assigning genes without overlaps
  * @return
  */
def getTopDomainWidth(x: Vector[LightBed], qtop: Double = 0.05,
  defaultScore: Double = -1.0): Double = {
  val widthsrt =
    listAllUniqueElement(x, defaultScore).map(y => y.score).sorted
  val i: Int = (widthsrt.length * (1 - qtop)).toInt
  widthsrt(i)
}

/** Get binary matrix for top domains.
  *
  * @param x
  *   matrix of regions, celltype by gene
  * @param width
  *   top domain definition
  * @return
  */
def getIsTopDomainMat(x: Vector[LightBed],
  width: Double): Vector[Vector[Boolean]] = {
  x.map(y => y.map(g => g.score >= width))
}

/** Get RTS (repressive tendency score) for the genes
  * Ref:
  * https://linkinghub.elsevier.com/retrieve/pii/S2405471220304191
  * @param x
  *   matrix of regions, celltype by gene domains
  * @param width
  *   top domain definition
  * @param defaultScore
  *   previously used for assigning the genes if no
  *   overlap found.
  * @return
  */
def getRTS(x: Vector[LightBed], width: Double,
  defaultScore: Double = -1.0): Vector[Double] = {
  val x0        = x.map(y => y.map(_.score))
  val nCellType = x.length
  val nGene     = x(0).length
  // get h_i,j
  val Bi =
    x0.map(y => mean(y.filter(x => x > defaultScore).toArray))
  val stdi =
    x0.map(y => stdev(y.filter(x => x > defaultScore).toArray))
  if (stdi.filter(x => x < 1e-10).nonEmpty) {
    println("Warning: some row(s) have extremely small std!")
    println("Remove those rows firstly.")
  } else {
    val meanBi   = mean(Bi.toArray)
    val meanStdi = mean(stdi.toArray)
    println(s"Across $nCellType celltypes: ")
    println(s"meanDomainWidth: $meanBi.")
    println(s"meanStd: $meanStdi.")
  }

  val h: Vector[Vector[Double]] = x0.zipWithIndex.map((y, i) => {
    y.map(b => (b - Bi(i)) / stdi(i))
  })
  // get v_i
  val v: Vector[Double] = Seq
    .range(0, nGene)
    .map(i => {
      h.map(hh => hh(i)).sum
    })
    .toVector

  // get v'_i
  val vmin                   = v.min
  val vmax                   = v.max
  val radius                 = vmax - vmin
  val vprime: Vector[Double] = v.map(x => (x - vmin) / radius)
  // get a_i
  val x1 = getIsTopDomainMat(x, width)
  val xi: Vector[Double] = Seq
    .range(0, nGene)
    .map(i => {
      x1.map(xx => ifelse(xx(i), 1.0, 0.0) + 1.0)
        .sum / (nCellType.toDouble + 1.0)
    })
    .toVector
  val a = vprime.zip(xi).map((v, x) => v * x)
  // get Ri
  val amin    = a.min
  val amax    = a.max
  val aRadius = amax - amin
  a.map(x => (x - amin) / aRadius)
}
