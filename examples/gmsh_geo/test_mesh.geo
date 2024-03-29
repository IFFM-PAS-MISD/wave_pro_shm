SetFactory("OpenCASCADE");
// Units [mm]

// plate
plate_thickness = 3;
L=40; // plate length
W=40; // plate width
plate_el_size = 3; //expected el size in mm

// delamination
// delamination semi-major and semi-minor axis
adelam=5; 
bdelam=5; 
// delamination centre coordinates
xdelam=L/4; 
ydelam=3*W/4; 
// delamination rotation
alpha=0*Pi; 
delam_thickness_pos = 1/3; // relative delamination position within plate thickness
below_delam_el = 1; // number of elements below delamination
above_delam_el = 2; // number of elements above delamination
delam_el_div = 25; // number of nodes on the outer edge of delamination (for Transfinite Curve)

// pzts
pzt_thickness = 0.5;
pzt_radius = 5;
xpzt1 = L/2;
ypzt1 = W/2;
xpzt2 = L/4;
ypzt2 = W/4;
pzt_el_div = 17; // number of nodes on the outer edge of pzt (for Transfinite Curve)
//

// piezoelectric transducers
// pzt1
pzt1_circle = newc; Circle(pzt1_circle) = {xpzt1, ypzt1, 0, pzt_radius, 0, 2*Pi};
Transfinite Curve{pzt1_circle} = pzt_el_div Using Progression 1;
pzt1_cl = newcl; Curve Loop(pzt1_cl) = {pzt1_circle};
pzt1_surf = news; Plane Surface(pzt1_surf) = {pzt1_cl}; 
// pzt2
pzt2_circle = newc; Circle(pzt2_circle) = {xpzt2, ypzt2, 0, pzt_radius, 0, 2*Pi};
Transfinite Curve{pzt2_circle} = pzt_el_div Using Progression 1;
pzt2_cl = newcl; Curve Loop(pzt2_cl) = {pzt2_circle};
pzt2_surf = news; Plane Surface(pzt2_surf) = {pzt2_cl}; 

// delamination
delam1_ellipse = newc;
Ellipse(delam1_ellipse) = {xdelam, ydelam, 0, adelam, bdelam, 0, 2*Pi};
//+
Rotate {{0, 0, 1}, {xdelam, ydelam, 0}, alpha} {
  Curve{delam1_ellipse}; 
}
Transfinite Curve{delam1_ellipse} = delam_el_div Using Progression 1;
delam1_cl = newcl; Curve Loop(delam1_cl) = {delam1_ellipse};
delam1_surf = news; Plane Surface(delam1_surf) = {delam1_cl}; 

// plate surf
p7 = newp; Point(p7) = {0,  0,  0,  1} ;
p8 = newp; Point(p8) = {L,  0,  0,  1} ;
p9 = newp; Point(p9) = {L,  W,  0,  1} ;
p10 = newp; Point(p10) = {0,  W,  0,  1} ;
l10 = newl; Line(l10) = {p7, p8};
l11 = newl; Line(l11) = {p8, p9};
l12 = newl; Line(l12) = {p9, p10};
l13 = newl; Line(l13) = {p10, p7};
Transfinite Curve{l10,l12} = L/plate_el_size Using Progression 1;
Transfinite Curve{l11,l13} = W/plate_el_size Using Progression 1;
plate_cl = newcl; Curve Loop(plate_cl) = {l10,l11,l12,l13};
plate_surf = news; Plane Surface(plate_surf) = {plate_cl,pzt1_cl,pzt2_cl,delam1_cl}; 

// Extrude pzts to inside of the plate considering delamination position
// pzt1
underpzt1_vol_below[] = Extrude {0, 0, delam_thickness_pos*plate_thickness} {
Surface{pzt1_surf}; Layers{below_delam_el}; Recombine;
};
underpzt1_vol_above[] = Extrude {0, 0, (1-delam_thickness_pos)*plate_thickness} {
Surface{underpzt1_vol_below[0]}; Layers{above_delam_el}; Recombine;
};
// pzt2
underpzt2_vol_below[] = Extrude {0, 0, delam_thickness_pos*plate_thickness} {
Surface{pzt2_surf}; Layers{below_delam_el}; Recombine;
};
underpzt2_vol_above[] = Extrude {0, 0, (1-delam_thickness_pos)*plate_thickness} {
Surface{underpzt2_vol_below[0]}; Layers{above_delam_el}; Recombine;
};
//Printf("underpzt1_vol_above[0] %g underpzt1_vol_above[1] %g underpzt1_vol_above[2] %g", underpzt1_vol_above[0], underpzt1_vol_above[1],underpzt1_vol_above[2]);

// Extrude delamination
// delamination 1
delam1_vol_below[] = Extrude {0, 0, delam_thickness_pos*plate_thickness} {
Surface{delam1_surf}; Layers{below_delam_el}; Recombine;
};
delam1_vol_above[] = Extrude {0, 0, (1-delam_thickness_pos)*plate_thickness} {
Surface{delam1_vol_below[0]}; Layers{above_delam_el}; Recombine;
};

// Extrude plate 
plate_vol_below[] = Extrude {0, 0, delam_thickness_pos*plate_thickness} {
Surface{plate_surf}; Layers{below_delam_el}; Recombine;
};
plate_vol_above[] = Extrude {0, 0, (1-delam_thickness_pos)*plate_thickness} {
Surface{plate_vol_below[0]}; Layers{above_delam_el}; Recombine;
};
//Printf("delam1_vol_below[0] %g delam1_vol_below[1] %g delam1_vol_below[2] %g", delam1_vol_below[0], delam1_vol_below[1],delam1_vol_below[2]);

//+ piezoelectric transducer volumes
// pzt1
pzt1_vol[] = Extrude {0, 0, pzt_thickness} {
  Surface{underpzt1_vol_above[0]}; 
Layers {1}; Recombine;
};
// pzt2
pzt2_vol[] = Extrude {0, 0, pzt_thickness} {
  Surface{underpzt2_vol_above[0]}; 
Layers {1}; Recombine;
};

// Additional points in surfaces for outputs
pzt1_center = newp; Point(pzt1_center) = {xpzt1,ypzt1,plate_thickness, 1} ;
Point{pzt1_center} In Surface {underpzt1_vol_above[0]};
pzt2_center = newp; Point(pzt2_center) = {xpzt2,ypzt2,plate_thickness, 1} ;
Point{pzt2_center} In Surface {underpzt2_vol_above[0]};

delam1_center = newp; Point(delam1_center) = {xdelam,ydelam,plate_thickness, 1} ;
Point{delam1_center} In Surface {delam1_vol_above[0]};

// Physical names
Physical Volume("plate") = {plate_vol_below[1],plate_vol_above[1],underpzt1_vol_below[1],underpzt1_vol_above[1],underpzt2_vol_below[1],underpzt2_vol_above[1]};
Physical Volume("pzt1") = {pzt1_vol[1]};
Physical Volume("pzt2") = {pzt2_vol[1]};
Physical Volume("delamination") = {delam1_vol_below[1],delam1_vol_above[1]};

Physical Surface("pzt1_electrode_plus") = {pzt1_vol[0]};
Physical Surface("pzt1_electrode_minus") = {underpzt1_vol_above[0]};
Physical Surface("pzt2_electrode_plus") = {pzt2_vol[0]};
Physical Surface("pzt2_electrode_minus") = {underpzt2_vol_above[0]};
Physical Surface("delamination_interface") = {delam1_vol_below[0]};

Physical Surface("plate_bottom_surf") = {plate_surf,pzt1_surf,pzt2_surf,delam1_surf};
Physical Surface("plate_top_surf") = {plate_vol_above[0],underpzt1_vol_above[0],underpzt2_vol_above[0],delam1_vol_above[0]};

Physical Curve("plate_left_top_edge") = {47};
Physical Curve("delamination_boundary") = {23};

Physical Point("pzt1_center") = {pzt1_center};
Physical Point("pzt2_center") = {pzt2_center};
Physical Point("delam1_center") = {delam1_center};
// Mesh parameters
// 1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms, 11: Quasi-structured Quad
Mesh.Algorithm = 6; // Frontal 6 default
// 1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT
Mesh.Algorithm3D = 1; // Delaunay 1 default
// 0: simple, 1: blossom, 2: simple full-quad, 3: blossom full-quad
Mesh.RecombinationAlgorithm = 1;// blossom 1 defauls
Mesh.CharacteristicLengthFactor = 1.00000;
Mesh.CharacteristicLengthMin = 0.10000;
//Mesh.CharacteristicLengthMax = 2.000000;
Mesh.CharacteristicLengthMax = 10.000000;
Mesh.RecombineAll = 1; // Apply recombination algorithm
Mesh.SubdivisionAlgorithm = 1; // 1: all quadrangles 2: all hexahedra, 3: barycentric
Mesh.Smoothing = 1;

Mesh.ElementOrder = 3;
Mesh 3;  // Generate 3D mesh
// Remove all duplicate mesh nodes in the current mesh:
Coherence Mesh;

// Mesh.SaveAll = 0; // check if this is the option for saving mesh - it must be 0 to output Physical groups//+

