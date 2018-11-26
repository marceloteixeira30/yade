/* lucScholtes2010 */

#pragma once

#include<pkg/common/ElastMat.hpp>
#include<pkg/common/Dispatching.hpp>
#include<pkg/common/NormShearPhys.hpp>
#include<pkg/dem/ScGeom.hpp>

/** This class holds information associated with each body state*/
class BPMpmState: public State {
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(BPMpmState,State,"JCFpm state information about each body.",
		((int,nbInitBonds,0,,"Number of initial bonds. [-]"))
		((int,nbBrokenBonds,0,,"Number of broken bonds. [-]"))
		((Real,damageIndex,0,,"Ratio of broken bonds over initial bonds. [-]"))
		,
		createIndex();
	);
	REGISTER_CLASS_INDEX(BPMpmState,State);
};
REGISTER_SERIALIZABLE(BPMpmState);

/** This class holds information associated with each body */
class BPMpmMat: public FrictMat {
  	public:
		virtual shared_ptr<State> newAssocState() const { return shared_ptr<State>(new BPMpmState); }
		virtual bool stateTypeOk(State* s) const { return (bool)dynamic_cast<BPMpmState*>(s); }
		
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(BPMpmMat,FrictMat,"Possibly jointed, cohesive frictional material, for use with other JCFpm classes",
		((int,type,0,,"If particles of two different types interact, it will be with friction only (no cohesion).[-]"))
		((Real,normalCohesion,0.,,"Defines the maximum admissible normal stress in traction in the matrix. [Pa]"))
		((Real,shearCohesion,0.,,"Defines the maximum admissible tangential stress in the matrix. [Pa]"))
		((Real,cohesiveYoung,0.,,"Defines the Young's modulus of the cohesive contact. [Pa]"))
		((Real,cohesivePoisson,0.,,"Defines the Poisson coefficient of the cohesive contact. [-]"))
		((Real,lambdaCohesion,0.,,"Parameter defining the radius of the cohesive contact. [-]"))
		((Real,betaCohesion,0.,,"Parameter to define the influence of the moments on the calculated maximum stress. [-]"))
                ((Real,residualFrictionAngle,-1.,,"Defines the residual friction angle (when contacts are not cohesive). residualFrictionAngle=frictionAngle if not specified. [degrees]"))
		,
		createIndex();
	);
	REGISTER_CLASS_INDEX(BPMpmMat,FrictMat);
};
REGISTER_SERIALIZABLE(BPMpmMat);

/** This class holds information associated with each interaction */
class BPMpmPhys: public NormShearPhys {
	public:
		virtual ~BPMpmPhys();

		YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(BPMpmPhys,NormShearPhys,"Representation of a single interaction of the JCFpm type, storage for relevant parameters",
			((Real,initD,0.,,"equilibrium distance for interacting particles. Computed as the interparticular distance at first contact detection."))
			((Real,beamNormalForce,0.,,"save the normal force for next step increment.[N]"))
			((Vector3r,beamShearForce,Vector3r::Zero(),,"save the shear force for next step increment.[N]"))
			((Real,beamNormalStiffness,0.,,"normal stiffness for the beam (cohesion between particles).[Pa]"))
			((Real,beamShearStiffness,0.,,"shear stiffness for the beam (cohesion between particles).[Pa]"))
			((Real,beamNormalCohesion,0.,,"Defines the maximum admissible normal stress in traction in the matrix. [Pa]"))
			((Real,beamShearCohesion,0.,,"Defines the maximum admissible tangential stress in the matrix for zero normal force. [Pa]"))
			((Real,beamSCoh,0.,,"Defines the maximum admissible tangential stress in the matrix for an applied normal force. [Pa]"))
			((Vector3r,beamMomentTwist,Vector3r(0,0,0),,"twist moment calculated for the beam.[Pa.m]"))
			((Vector3r,beamMomentBending,Vector3r(0,0,0),,"bending moment calculated for the beam.[Pa.m]"))
			((bool,isBroken,false,,"flag for broken interactions"))
			((bool,isCohesive,false,,"If false, particles interact in a frictional way. If true, particles are bonded regarding the given :yref:`cohesion<JCFpmMat.cohesion>` and :yref:`tensile strength<JCFpmMat.tensileStrength>` (or their jointed variants)."))
			((Real,tanFrictionAngle,0.,,"tangent of Coulomb friction angle for this interaction (auto. computed). [-]"))
			((Real,tanResidualFrictionAngle,0.,,"tangent of Coulomb residual friction angle for this interaction (friction contact). [-]"))
			((Real,beamRadius,0.,,"beamRadius=lambda*Rmin. [m]"))
			((Real,beamArea,0.,,"beamArea=pi*beamRadius^2. [m^2]"))
			((Real,beamMomInertia,0.,,"beamMomInertia=(1/4)*pi*R^4. [m^4]"))
			((Real,beamPolarMomInertia,0.,,"beamPolarMomInertia=(1/2)*pi*R^4. [m^4]"))
			((Real,previousDisplacement,0.,,"displacement between two particles from the previous step. [m]"))
			((bool,breakOccurred,0,,"Flag used to trigger retriangulation as soon as a cohesive bond breaks in FlowEngine (for DFNFlow use only)"))
			((Real,beamBeta,0.,,"Parameter to define the influence of the moments on the calculated maximum stress. [-]"))
			,
			createIndex();
			,
		);
		DECLARE_LOGGER;
		REGISTER_CLASS_INDEX(BPMpmPhys,NormShearPhys);
};
REGISTER_SERIALIZABLE(BPMpmPhys);

/** 2d functor creating InteractionPhysics (Ip2) taking JCFpmMat and JCFpmMat of 2 bodies, returning type JCFpmPhys */
class Ip2_BPMpmMat_BPMpmMat_BPMpmPhys: public IPhysFunctor{
	public:
		virtual void go(const shared_ptr<Material>& pp1, const shared_ptr<Material>& pp2, const shared_ptr<Interaction>& interaction);
		
		FUNCTOR2D(BPMpmMat,BPMpmMat);
		DECLARE_LOGGER;
		YADE_CLASS_BASE_DOC_ATTRS(Ip2_BPMpmMat_BPMpmMat_BPMpmPhys,IPhysFunctor,"Converts 2 :yref:`JCFpmMat` instances to one :yref:`JCFpmPhys` instance, with corresponding parameters. See :yref:`JCFpmMat` and [Duriez2016]_ for details",                   
			((int,cohesiveTresholdIteration,1,,"should new contacts be cohesive? If strictly negativ, they will in any case. If positiv, they will before this iter, they won't afterward."))
		);
		
};
REGISTER_SERIALIZABLE(Ip2_BPMpmMat_BPMpmMat_BPMpmPhys);

/** 2d functor creating the interaction law (Law2) based on SphereContactGeometry (ScGeom) and JCFpmPhys of 2 bodies, returning type JointedCohesiveFrictionalPM */
class Law2_ScGeom_BPMpmPhys_BondedContactM: public LawFunctor{
	public:
		virtual bool go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I);
		FUNCTOR2D(ScGeom,BPMpmPhys);

		YADE_CLASS_BASE_DOC_ATTRS(Law2_ScGeom_BPMpmPhys_BondedContactM,LawFunctor,"Interaction law for cohesive frictional material, e.g. rock, possibly presenting joint surfaces, that can be mechanically described with a smooth contact logic [Ivars2011]_ (implemented in Yade in [Scholtes2012]_). See examples/jointedCohesiveFrictionalPM for script examples. Joint surface definitions (through stl meshes or direct definition with gts module) are illustrated there.",
			((bool,neverErase,false,,"Keep interactions even if particles go away from each other (only in case another constitutive law is in the scene"))
			((bool,cracksFileExist,false,,"if true (and if :yref:`recordCracks<Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM.recordCracks>`), data are appended to an existing 'cracksKey' text file; otherwise its content is reset."))
			((string,Key,"",,"string specifying the name of saved file 'cracks___.txt', when :yref:`recordCracks<Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM.recordCracks>` is true."))
			((bool,recordCracks,false,,"if true, data about interactions that lose their cohesive feature are stored in the text file cracksKey.txt (see :yref:`Key<Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM.Key>` and :yref:`cracksFileExist<Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM.cracksFileExist>`). It contains 9 columns: the break iteration, the 3 coordinates of the contact point, the type (1 means shear break, while 0 corresponds to tensile break), the ''cross section'' (mean radius of the 2 spheres) and the 3 coordinates of the contact normal."))
			((int,nbTensCracks,0,,"number of tensile microcracks."))
			((int,nbShearCracks,0,,"number of shear microcracks."))
		);
		DECLARE_LOGGER;	
};
REGISTER_SERIALIZABLE(Law2_ScGeom_BPMpmPhys_BondedContactM);
