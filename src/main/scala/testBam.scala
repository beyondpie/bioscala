import os._
import htsjdk.* 
import htsjdk.samtools.DefaultSAMRecordFactory
import htsjdk.samtools.SAMFileWriter
import htsjdk.samtools.SAMFileWriterFactory
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamInputResource
import htsjdk.samtools.SamReader
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.ValidationStringency
import htsjdk.samtools.seekablestream.SeekableStream

import java.io.File
import java.io.IOException
import java.net.MalformedURLException
import java.net.URL



def getnread(fnm: String): Int = {
  val fnm = os.pwd / "src" / "test" / "resource" / "H3K27ac-Male.srt.bam"
  // header does not count
  var n = -1
  val reader =
    SamReaderFactory.makeDefault().open(new File(fnm.toString))

  val it = reader.iterator()

  while (it.hasNext) {
    it.next
    n = n + 1
  }
  reader.close()
  n
}

