// run with
// openscad -o robin.stl merge.scad

union() {
  import("fuselage.stl");
  import("pylon.stl");
}
