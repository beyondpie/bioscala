package ZSTCompress

import SZUtils.TaskElement

class gz2zst(
  val gzf: String, val zstf: String, val logfnm: String,
  val flagfnm: String, val skip: Boolean = true,
  val zstd: String = "/home/szu/miniforge3/bin/zstd",
  val T: Int = 10
) extends TaskElement {

  def runCore(): Unit = {
    if (!os.exists(os.Path(zstf))) {
      println(s"gz2zst for ${zstf}")
      os.proc("gzip", "-d", "-c", gzf)
        .pipeTo(os.proc(zstd, s"-T${T}", "-19", "-f", "-o", zstf))
        .call(stdout = os.Path(logfnm), stderr = os.Path(logfnm),
            check = true)
      println(s"gz2zst for ${zstf} done.")
    } else {
      println(s"${zstf} exists, and skip it.")
    }
  }
}
