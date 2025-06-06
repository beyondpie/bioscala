// set to Debug for compilation details (Info is default)
logLevel := Level.Info

// main
ThisBuild / version := "0.0.33"
ThisBuild / organization := "io.github.beyondpie"
ThisBuild / organizationName := "zulab"
ThisBuild / scalaVersion := "3.7.0"
name  := "bioscala"

// for igv
resolvers += "Bioviz" at "https://nexus.bioviz.org/repository/maven-releases/"

// dependencies
libraryDependencies += "org.scalatest" %% "scalatest" % "3.2.19" % "test"
libraryDependencies += "com.lihaoyi" %% "os-lib" % "0.11.4"
libraryDependencies += "org.scala-lang.modules" %% "scala-parallel-collections" % "1.2.0"
// libraryDependencies += "org.plotly-scala" % "plotly-render_2.13" % "0.8.5"
libraryDependencies += "com.github.samtools" % "htsjdk" % "4.2.0"
// libraryDependencies += "com.lihaoyi" % "ammonite" % "3.0.1" cross CrossVersion.full
libraryDependencies += "commons-io"           % "commons-io"  % "2.19.0"
libraryDependencies += "org.jsoup"            % "jsoup"       % "1.20.1"
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
libraryDependencies += "org.slf4j" % "slf4j-simple" % "2.0.17"
libraryDependencies += "org.slf4j" % "slf4j-api"    % "2.0.17"
libraryDependencies += "org.broad.igv" % "bigwig" % "3.0.0"

libraryDependencies ++= Seq(
 "dev.optics" %% "monocle-core"  % "3.3.0",
 "dev.optics" %% "monocle-macro" % "3.3.0",
)

// scalaTest settings
logBuffered in Test := false

scalacOptions ++= Seq(
    "-encoding",
    "utf8",
    "-feature",
    "-language:implicitConversions",
    "-language:existentials",
    "-experimental",
    "-unchecked",
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

ThisBuild / organizationHomepage := Some(url("https://beyondpie.github.io/"))

ThisBuild / scmInfo := Some(
  ScmInfo(
    url("https://github.com/beyondpie/bioscala"),
    "scm:git@github.com:beyondpie/bioscala.git"
  )
)

ThisBuild / versionScheme := Some("early-semver")

ThisBuild / developers := List(
  Developer(
    id = "io.github.beyondpie",
    name = "Songpeng Zu",
    email = "panshi.pie@email",
    url = url("https://beyondpie.github.io/")
  )
)

ThisBuild / description := "Explore the usage of scala on biology"
ThisBuild / licenses := List(
  "GPLv3" -> new URL("https://www.gnu.org/licenses/gpl-3.0.en.html")
)
ThisBuild / homepage := Some(url("https://github.com/beyondpie/bioscala"))

// Remove all additional repository other than Maven Central from POM
ThisBuild / pomIncludeRepository := { _ => false }
ThisBuild / publishMavenStyle := true


// new setting for the Central Portal
ThisBuild / publishTo := {
  val centralSnapshots = "https://central.sonatype.com/repository/maven-snapshots/"
  if (isSnapshot.value) Some("central-snapshots" at centralSnapshots)
  else localStaging.value
}
