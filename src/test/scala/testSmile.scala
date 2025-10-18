import smile.math.special.*
import smile.*
import smile.clustering.*
import smile.clustering.linkage.*

object TestSmile {
  Erf.erf(1.0)
  val projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
  val x = read.csv(s"${projd}/src/test/resource/rem.txt",
     header=false, delimiter=" ").toArray()
  val clusters = hclust(x, "complete")
  val y = clusters.partition(6)
}
