// set to Debug for compilation details (Info is default)
logLevel := Level.Info

// main
ThisBuild / version := "0.0.99"
ThisBuild / organization := "com.szu"
ThisBuild / scalaVersion := "3.6.3"
name  := "bioscala"

// for igv
resolvers += "Bioviz" at "https://nexus.bioviz.org/repository/maven-releases/"

// dependencies
libraryDependencies += "com.lihaoyi" %% "os-lib" % "0.11.3"
libraryDependencies += "org.scala-lang.modules" %% "scala-parallel-collections" % "1.2.0"
// libraryDependencies += "org.plotly-scala" % "plotly-render_2.13" % "0.8.5"
libraryDependencies += "com.github.samtools" % "htsjdk" % "4.1.1"
libraryDependencies += "com.lihaoyi" % "ammonite" % "3.0.1" cross CrossVersion.full
libraryDependencies += "commons-io"           % "commons-io"  % "2.17.0"
libraryDependencies += "org.jsoup"            % "jsoup"       % "1.18.3"
libraryDependencies += "com.github.haifengl" %% "smile-scala" % "4.2.0"
libraryDependencies ++= Seq(
    "org.bytedeco" % "javacpp" % "1.5.11"
      classifier "macosx-arm64" classifier "macosx-x86_64"
      classifier "windows-x86_64" classifier "linux-x86_64",
    "org.bytedeco" % "openblas" % "0.3.28-1.5.11"
      classifier "macosx-arm64" classifier "macosx-x86_64"
      classifier "windows-x86_64" classifier "linux-x86_64",
    "org.bytedeco" % "arpack-ng" % "3.9.1-1.5.11"
      classifier "macosx-x86_64" classifier "windows-x86_64"
      classifier "linux-x86_64"
)
libraryDependencies += "org.slf4j" % "slf4j-simple" % "2.0.16"
libraryDependencies += "org.slf4j" % "slf4j-api"    % "2.0.16"
libraryDependencies += "org.broad.igv" % "bigwig" % "3.0.0"

libraryDependencies ++= Seq(
 "dev.optics" %% "monocle-core"  % "3.1.0",
 "dev.optics" %% "monocle-macro" % "3.1.0",
)

scalacOptions ++= Seq(
    "-encoding",
    "utf8",
    "-feature",
    "-language:implicitConversions",
    "-language:existentials",
    "-unchecked",
    "-Werror",
    "-explain-types",
    "-explain",
    "-deprecation"
)

// * to generate assembled jar
// This will slow down REPL start a lot.

// Compile / packageBin / mappings := {
//   val dependencyJars = (Compile / dependencyClasspath).value.files
//   val compiledClasses = (Compile / packageBin / mappings).value
//   compiledClasses ++ dependencyJars.map(jar => jar -> s"lib/${jar.getName}")
// }
// assembly / assemblyMergeStrategy := {
//   case PathList("META-INF", xs @ _*) => MergeStrategy.discard
//   case x => MergeStrategy.first
// }
