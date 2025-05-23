package bioscala.TRIAGE

import Genome.MouseGenome.chr2size
import SZUtils.ifelse
import bioscala.LightCoord.Coord
import bioscala.LightCoord.Coords
import bioscala.LightCoord.GFF3
import bioscala.LightCoord.GeneBody3to5
import bioscala.LightCoord.GeneBody5to3
import bioscala.LightCoord.LightBed
import bioscala.LightCoord.LightBedElement
import bioscala.LightCoord.aroundGene
import bioscala.LightCoord.coordOrd
import bioscala.LightCoord.findOvlpOneChrSorted
import bioscala.LightCoord.genomeCoordIgnoreStrandOrd
import bioscala.LightCoord.getCenter
import smile.math.MathEx.mean
import smile.math.MathEx.stdev

import math.*

/**
 * Group Bed elements into chromosome majored, and sort the
 * coords within each chromosome.
 *
 * @param d vector of LightBedElements
 * @return
 */
def groupSortLightBed(d: LightBed): Map[String, LightBed] = {
  d.groupBy(_.g.chrom)
    .map((chr, v) => (chr, v.sortBy(_.g.coord)(using coordOrd)))
}

/**
 * Assign Broadest domain to gene bodies.
 * Ref: TRIAGE, https://linkinghub.elsevier.com/retrieve/pii/S2405471220304191
 *
 * For each gene body, a LightBedElement will be got. If no broad domain found, will
 * assign one with genome coords as extended gene body,  name as "None",
 * score (domain width) as 0.0.
 *
 * 1. Extend the gene body based on upStream to TSS, and downStream to TES
 * as proximal region, see [[aroundGene]].
 * 2. Map the center / any position of broad domain to the proximal region.
 * 3 .Only broadest domain used for each gene.
 *
 * @param x genes
 * @param d broad domain
 * @param useDomainCenter if only domain center considered for overlapping
 *  In TRIGAGE, it seems that they firslty use domain center for each gene's overlapping,
 *  Then later they consider all the domains having at least one gene mapped above, and
 *  check if the domains also overlap with other genes though no center overlapped.
 * @return
 */
def assignBroadestDomain2Gene(genes: GFF3, d: LightBed): LightBed = {
  val dsrt = groupSortLightBed(d)
  genes.zipWithIndex
    .groupBy(_._1.g.chrom)
    .flatMap((chrom, gs) => {
      if (!dsrt.contains(chrom)) {
        gs.map(t => ((g = t._1.g, name = "None", score = 0.0), t._2))
      } else {
        // prepare inputs for findOvlpOneChrSorted
        val gsrt: Vector[(Coord, Int)] = gs
          .map(x => (x._1.g.coord, x._2))
          .sortBy(_._1)(using coordOrd)

        // run it to assign some gene bodies broad domains.
        val gIndex2d: Map[Int, LightBedElement] =
          findOvlpOneChrSorted(query = gsrt.map(_._1),
              subject = dsrt(chrom).map(_.g.coord))
            .map(x => {
              val dbrdst: LightBedElement = Seq
                .range(x._2.startFrom, x._2.endTo)
                .map(i => dsrt(chrom)(i))
                .maxBy(_.score) // only broadest domain kept
              (gsrt(x._1)._2, dbrdst)
            })
            .toMap
        // assign each gene body a broad domain
        gsrt
          .map(_._2)
          .map(i => {
            // if not matched, use "None" with gene's coord
            val r = gIndex2d.getOrElse(i,
                (g = genes(i).g, name = "None", score = 0.0))
            (r, i)
          })
      }
    })
    .toVector     // flatMap to Vector
    .sortBy(_._2) // keep the input order by sorting the index
    .map(_._1)    // remove index
}

/**
 * Get all the unique LightBedElements.
 * The elements will be ordered by genomeCoordIgnoreStrandOrd
 * @param x matrix of LightBedElement
 * @param minScore filter element by minScore before listing all of them
 * @return
 */
def listAllUniqueElement(x: Vector[LightBed],
  minScore: Double = 0.0): LightBed = {
  x.flatMap(y => y.filter(g => g.score >= minScore))
    .distinct
    .sortBy(_.g)(using genomeCoordIgnoreStrandOrd)
}

/**
 * Get the top domain width from the mapped domains
 *
 * Ref: https://linkinghub.elsevier.com/retrieve/pii/S2405471220304191
 *
 * @param x assigned domains organized as subclass by genes
 * @param qtop top ratio, which should be determined by the coverages of
 *  variably expressed TFs (VETFs) or other genes using the top domains.
 *  Follow the reference, we use 0.05 as default.
 * @return
 */
def getTopDomainWidth(x: Vector[LightBed],
  qtop: Double = 0.05): Double = {
  val widthsrt = listAllUniqueElement(x).map(y => y.score).sorted
  val i: Int   = widthsrt.length * (1 - qtop).toInt
  widthsrt(i)
}

/**
 * Get binary matrix for top domains.
 *
 * @param x matrix of regions
 * @param width top domain definition
 * @return
 */
def getIsTopDomainMat(x: Vector[LightBed],
  width: Double): Vector[Vector[Boolean]] = {
  x.map(y => y.map(g => g.score >= width))
}

/**
 * Get RTS (repressive tendency score) for the genes
 * Ref: https://linkinghub.elsevier.com/retrieve/pii/S2405471220304191
 * @param x matrix of regions, celltype by mapped broad domains
 * @param width top domain definition
 * @return
 */
def getRTS(x: Vector[LightBed], width: Double): Vector[Double] = {
  val x0        = x.map(y => y.map(_.score))
  val nCellType = x.length
  val nGene     = x(0).length
  // get h_i,j
  val Bi   = x0.map(y => mean(y.toArray))
  val stdi = x0.map(y => stdev(y.toArray))
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
        .sum / nCellType.toDouble
    })
    .toVector
  val a = vprime.zip(xi).map((v, x) => v * x)
  // get Ri
  val amin    = a.min
  val amax    = a.max
  val aRadius = amax - amin
  a.map(x => (x - amin) / aRadius)
}

