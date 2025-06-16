package bioscala.LightCoord.BedGraph

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

// TODO: make it as a formal test suite.
class LightBedGraphTest {
  // TODO: replace with project related automatical path..
  val projd =
    "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val bg2bw =
    "/home/szu/miniforge3/envs/scala/bin/bedGraphToBigWig"

  val workd = os.Path(projd) / "bioscala"
  val resourced = workd / "src"/ "test" / "resources"
  val bedgraphfnm = resourced / "test.bedgraph"
  val chromsizefnm = resourced / "test.chromSize.txt"
  val bwfnm = resourced / "test.bigwig"

  // * test transform to bigwig from bedgraph
  LightBedGraph.toBigWig(bg2bw, chromsizefnm,
    bedgraphfnm, bwfnm)

  // * test reader
  val r = getBigWigReader(bwfnm)
  val x = (chrom = "chr1",
    coord = (startFrom = 100, endTo = 101), strand = ".")
  val y = getWigValue(r, x,
    emptyValue = -1.0, contained = false)

}
