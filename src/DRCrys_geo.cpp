//==========================================================================
//  AIDA Detector description implementation 
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Frank
//
//==========================================================================
//
// Specialized generic detector constructor
// 
//==========================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "DD4hep/OpticalSurfaces.h"


using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {
  std::cout<<"Creating DRCrys "<<std::endl;
// static means only accessible in this file
//// Function ot create a detector.
//description: a reference to a detector object
// e: an xml node, which likely ocntains configuration data for the detector to be created.
// sens: SensitiveDetector obkect, likely used to specify how the detector should respond to particles

  static double tol = 0.001; //tolerance value, likely used for some comparison

  // material to underly it all
  Material      air       = description.air();
// function retrieves air material from detector object 

  // get stuff from xml file
  xml_det_t     x_det     = e;
  Layering      layering (e);
// gets data from xml file to define the crystal : check this out

  // for volume tags in detector
  int           det_id    = x_det.id();
  std::cout<<" det_id is "<<det_id<<std::endl;
  string        det_name  = x_det.nameStr();
  std::cout<<" det_name is .. "<<det_name<<std::endl;
//parses xml node, extracts detector ID and name
//

  // pointer to finding dimensins text in xml file
  // look in DDCore/include/Parsers/detail/Dimension.h for arguments
  xml_comp_t    x_towers  = x_det.staves(); //Retrieves a sub-component from the detector element x_det called "staves" - individual segments/posts?
  xml_comp_t    x_dim     = x_det.dimensions();// Retrieves the 'dimensions' subcomponent from detectorelement, contains info about size and placement?



  double hwidth   = x_dim.width()/2.;//calculates half-width based on width retrieved from x_dim
  double hzmax = x_dim.z_length()/2.;//half-length of det along z-axis
  int Ncount = x_dim.repeat();//Extracts # of repetitions of the detector element, maybe number of staves? 
  double zoffset = x_dim.z1();//z-offset from dimensions, used to position detector relative to a reference plane


  std::cout<<"half width zmax are "<<hwidth<<" "<<hzmax<<std::endl;
  std::cout<<" array size is "<<Ncount<<std::endl;
  std::cout<<" z offset is "<<zoffset<<std::endl;



  double agap = x_dim.gap();//gao size between elements from x_dim, likely used when placing multiple components 
  std::cout<<" gap between array elements is "<<agap<<std::endl;



  OpticalSurfaceManager surfMgr = description.surfaceManager();//Retrieves optical surface manager from the detector description, likely used to manage reflective of refractive properties for light interactions within the detector
  OpticalSurface cryS  = surfMgr.opticalSurface("/world/"+det_name+"#mirrorSurface");//Retrieves an optical surface names 'mirrorSurface' from the surface manager, used to model how light reflects within detector?



  
  // detector element for entire detector.  
  DetElement    sdet      (det_name,det_id);//created detector element object sdet using previously extracted detector name and ID
  Volume        motherVol = description.pickMotherVolume(sdet);//Retrieves the mother volume, likely existing volume within the detector setup


  //PolyhedraRegular hedra  (nphi,inner_r,outer_r+tol*2e0,zmaxt);
  //set containment area for whole calorimeter
  Box abox  ((2*Ncount+1)*(hwidth+agap+tol),(2*Ncount+1)*(hwidth+agap+tol), hzmax+tol);// Defines a box shape for hte detector envelope based on the number of repetitions, half-width, gap size, tolerance, and half-length
  Volume        envelope  (det_name,abox,air);//Creates a volume object named 'envelope' using defined box shape and air material
  Position a_pos(0.,0.,(hzmax+zoffset));//defines a position for the detector envelope at origin along x and y and offset along zaxis by half-length and z-offset values
  PlacedVolume  env_phv   = motherVol.placeVolume(envelope,a_pos);//places the envelope volume within the mother volume at the specified position a_pos. This creates a physical instance of the detector element in simulation

  env_phv.addPhysVolID("system",det_id);//assigns a physical volume id "system" along the detector ID to the placed volume, maybe used ot identify the detector element for tracking particles or assigning properties?

  sdet.setPlacement(env_phv);  // associate the placed volume to the detector element - creates a link between physical volume and logical det element
  sens.setType("calorimeter");//sets type of sensitive detector 'sens' to calorimeter - most likely means that detector measures energy deposited by particles



  // create towers to put into the calorimeter

  // a tower may have several different patterns that repeat.
  //   for example, there may be 10 layers with one thickness of Pb and scint and 20 with anothe r set of thicknesses.
  //   each of these repeating things is a "layer". (so in this example, two "layers")
  //within a layer is a slice os the Pb and scint are slices
  //the assembled tower is a Stave


    // tower envelope
  dd4hep::Box towertrap(hwidth+tol,hwidth+tol,hzmax+tol);// Defines a box shape for a single tower elemtn based on half width and half length
  dd4hep::Volume towerVol( "tower", towertrap, air);// Sets volume object using defined box shape and air material
  std::cout<<"   tower visstr is "<<x_towers.visStr()<<std::endl;//prints visualization string associated with the detector towers from the XML configuration. likely used to define how the detector element is displayed in the simulation software
  towerVol.setVisAttributes(description, x_towers.visStr());//sets visualization attributes for the tower volume based on the visualization string retrieved from the XML configuration
  towerVol.setSensitiveDetector(sens);//assigns the same sensitive detector 'sens' to the tower volume, indicating that individual towers are also sensitive to particle interactions 

  int itower=0; 
  string t_name1 = _toString(itower,"towerp%d") ;//sets a string t_name1 representing itower?
  DetElement tower_det(t_name1,det_id);  // detector element for a tower





  //  SkinSurface haha = SkinSurface(description,sdet, "HallCrys", cryS, towerVol);
  //  haha.isValid();



    // Loop over the sets of layer elements in the detector.
  double z_bottoml  = -hzmax;// set bottom of detector to -hzmax value
  int l_num = 1;
  for(xml_coll_t li(x_det,_U(layer)); li; ++li)  { // iterating over layers within detector element x_det
    std::cout<<"DRCrys layer (layers contain slices of material)"<<l_num<<std::endl;//indicates the start of processing a layer element
    xml_comp_t x_layer = li;// extracts current layer element from loop iteration 
    int repeat = x_layer.repeat(); //retrieves the number of times the current layer element should be repeated over the detector
      // Loop over number of repeats for this layer.
    for (int j=0; j<repeat; j++)    {//iterates through  # of reps
      std::cout<<"DRCrys layer "<<li<<" repeat "<<j<<std::endl;
      string l_name = _toString(l_num,"layer%d");//unique name for current layer instance
      double l_hzthick = layering.layer(l_num-1)->thickness()/2.;  // Layer's thickness. Ethan: should this not be half-thickness?
      std::cout<<"half  thickness is "<<l_hzthick<<std::endl;//okay yeah it is half thickness

	// find top and bottom lengths at this position and center
        // relative to tower bottom
      double z_topl=z_bottoml + 2.*l_hzthick;
      double z_midl=z_bottoml + l_hzthick;
      Position   l_pos(0.,0.,z_midl);      // Position of the layer. (along z axis)
      std::cout<<" placed at z of "<<z_midl<<std::endl;

      dd4hep::Box l_box(hwidth,hwidth,l_hzthick); // creates a box shape for current layer volume
      dd4hep::Volume     l_vol(l_name,l_box,air;//creates a volume object for current layer using defined box shape and air material - likely space occupied by the layer within the tower
      DetElement layer(tower_det, l_name, det_id);//creates a detector element object for crrent layer, associating it with tower detector element tower_det and using detector ID

        // Loop over the sublayers or slices for this layer.
      int s_num = 1;

      double z_bottoms2=-l_hzthick;  
      for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
	std::cout<<" with slice "<<s_num<<std::endl;
	xml_comp_t x_slice = si;
	string     s_name  = _toString(s_num,"slice%d");
	double     s_hzthick = x_slice.thickness()/2.;
	std::cout<<" with half  thickness "<<s_hzthick<<" and material "<<x_slice.materialStr()<<std::endl;

	      // this is relative to tower bottom, not layer bottom

	double z_mids2 = z_bottoms2+s_hzthick;
	      

	Position   s_pos(0.,0.,z_mids2);      // Position of the layer.
	std::cout<<" placed at "<<z_mids2<<std::endl;
	dd4hep::Box s_box(hwidth,hwidth,s_hzthick);


	dd4hep::Volume     s_vol(s_name,s_box,description.material(x_slice.materialStr()));//creates a volume object for the current slice using the decined box shape and the material manager of the detector description 
	DetElement slice(layer,s_name,det_id);//creates a detector element object for the current slice, associating it with the parent layer element 'layer' and using detector ID det_id

	if ( x_slice.isSensitive() ) { //if current slice element is flagged as sensitive in XML config...
		s_vol.setSensitiveDetector(sens);//... the same senitive detector 'sens' is assigned to its volume, indicating that this material component will also be sensitive to particle interactions 
          }
	std::cout<<"          slice visstr is "<<x_slice.visStr()<<std::endl;//retrieves and prints out visualization string for the slice
	slice.setAttributes(description,s_vol,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());//sets attributes for the slice detector element, including the region it covers within the detector, its limits/boundaries, and vis style 

          // Slice placement.
	PlacedVolume slice_phv = l_vol.placeVolume(s_vol,s_pos);//places the slice volume 's_vol' within the parent layer volyume 'l_vol' at the calculated position 's_pos', creating a physical instance of slice material
	slice_phv.addPhysVolID("slice", s_num);//assigns a physical volume ID along with the slice number to the placed volume 
	slice.setPlacement(slice_phv);//associates the placed volume iwth the slice detector element 
          // Increment Z position of slice.

	z_bottoms2 += 2.*s_hzthick;//updates the bottom position for the next slice iwthin the layer based on the current slice's thickness

          // Increment slice number.
	      ++s_num;
      }//end of loop over sublayers/slices for each layer


      // place the layer
        // Set region, limitset, and vis of layer.
      std::cout<<" layer visstr is "<<x_layer.visStr()<<std::endl;
      layer.setAttributes(description,l_vol,x_layer.regionStr(),x_layer.limitsStr(),x_layer.visStr());

      PlacedVolume layer_phv = towerVol.placeVolume(l_vol,l_pos);
      layer_phv.addPhysVolID("layer", l_num);
      layer.setPlacement(layer_phv);
        // Increment to next layer Z position.


      z_bottoml=z_bottoml+2.*l_hzthick;

      ++l_num;
    }
  }
      
  // now that you put the layers and slices into the tower, place it  

    // now do the placement




  int towernum=-1;
  for (int ijk1=-Ncount; ijk1<Ncount+1; ijk1++) {//iterate through number of repetitions from the detector config
    for (int ijk2=-Ncount; ijk2<Ncount+1; ijk2++) {//iterates over all possible combinations of indices, creating a grid of positions
      double mod_x_off = (ijk1)*2*(hwidth+tol+agap);//x offset             
      double mod_y_off = (ijk2)*2*(hwidth+tol+agap);//y offset
      std::cout<<"placing crystal at ("<<mod_x_off<<","<<mod_y_off<<")"<<std::endl;


      Transform3D tr(RotationZYX(0.,0.,0.),Position(mod_x_off,mod_y_off,0.));//defines transformation object tr that specifies position (with calculated offsets) and 0 rotation
      PlacedVolume pv = envelope.placeVolume(towerVol,tr);//places the tower volume towerVol within the envelope volume, creating a physical tower inside the detector



      pv.addPhysVolID("system",det_id);//assigns the same system ID as the overall detector to the placed tower volume
      pv.addPhysVolID("ix",ijk1);//assigns physical volume ID ijk1 as an identifier of the tower's position of x value
      pv.addPhysVolID("iy",ijk2);//for y-coordinate


      towernum +=1;//keep track of number of towers placed
      std::cout<<"placing tower "<<towernum<<std::endl;
      string t_name2 = _toString(towernum,"0%d");//generates a name for current tower element
      DetElement sd = tower_det.clone(t_name2,det_id);//creates a new detector element sd by cloning the parent tower_det and assigning new name and det ID



      sd.setPlacement(pv);//sets placement ??

      string tt_name = _toString(towernum,"HallCrys%d");//generates new name 
      BorderSurface haha = BorderSurface(description,sdet, tt_name, cryS, pv,env_phv); //creates border surface object (haha) using the detector description, sdet, the name, the optical surface cryS, the placed tower volume pv, and placed envelope volume env_phv
//relates how light interacts with the edges or boundaries of the tower
      haha.isValid();//checks the validity of the created border surface object




      //      sdet.add(sd);



  }//end of iterating over all positions ijk2

  }//end of iterating over all positions ijk1

  // Set envelope volume attributes.
      std::cout<<" envelope visstr is "<<x_det.visStr()<<std::endl;
    envelope.setAttributes(description,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());




  std::cout<<"exiting DRCrys creator"<<std::endl;

  return sdet;
}

DECLARE_DETELEMENT(DRCrys,create_detector)


