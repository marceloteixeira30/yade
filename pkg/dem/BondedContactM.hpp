/* MarceloDeSimone2018 */

#pragma once

#include<pkg/common/ElastMat.hpp>
#include<pkg/common/Dispatching.hpp>
#include<pkg/common/NormShearPhys.hpp>
#include<pkg/dem/ScGeom.hpp>

/** This class holds information associated with each body state*/
class BPMState: public State {
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(BPMState,State,"BPM state information about each body.",
		((int,nbInitBonds,0,,"Number of initial bonds. [-]"))
		((int,nbBrokenBonds,0,,"Number of broken bonds. [-]"))
		((Real,damageIndex,0,,"Ratio of broken bonds over initial bonds. [-]"))
		,
		createIndex();
	);
	REGISTER_CLASS_INDEX(BPMState,State);
};
REGISTER_SERIALIZABLE(BPMState);
JointedCohesiveFrictionalPM
/** This class holds information associated with each body */
class BPMMat: public FrictMat {
  	public:
		virtual shared_ptr<State> newAssocState() const { return shared_ptr<State>(new BPMState); }
		virtual bool stateTypeOk(State* s) const { return (bool)dynamic_cast<BPMState*>(s); }
		
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(BPMMat,FrictMat,"Cohesive frictional material, for use with other BPM classes",
		((int,type,0,,"If particles of two different types interact, it will be with friction only (no cohesion).[-]"))
		((Real,residualFrictionAngle,-1.,,"Defines the residual friction angle (when contacts are not cohesive). residualFrictionAngle=frictionAngle if not specified. [degrees]"))
		((Real,normalCohesion,0.,,"Defines the maximum admissible normal stress in traction in the matrix. [Pa]"))
		((Real,shearCohesion,0.,,"Defines the maximum admissible tangential stress, for Fn=0, in the matrix. [Pa]"))
		((Real,cohesiveYoung,0.,,"Defines the Young's modulus of the cohesive contact. [Pa]"))
		((Real,cohesivePoisson,0.,,"Defines the Poisson coefficient of the cohesive contact. [-]"))
		((Real,lambda,0.,,"Parameter defining the radius of the cohesive contact. [-]"))
		((Real,beta,0.,,"Parameter to define the influence of the moments on the calculated maximum stress. [-]"))
		,
		createIndex();
	);
	REGISTER_CLASS_INDEX(BPMMat,FrictMat);
};
REGISTER_SERIALIZABLE(BPMMat);

/** This class holds information associated with each interaction */
class BPMPhys: public NormShearPhys {
	public:
		virtual ~BPMPhys();

		YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(BPMPhys,NormShearPhys,"Representation of a single interaction of the BPM type, storage for relevant parameters",
			((Real,beamNormalForce,0.,,"save the normal force for next step increment.[N]"))
			((Real,beamShearForce,0.,,"save the shear force for next step increment.[N]"))
			((Real,beamNormalStiffness,0.,,"normal stiffness for the beam (cohesion between particles).[Pa]"))
			((Real,beamShearStiffness,0.,,"shear stiffness for the beam (cohesion between particles).[Pa]"))
			((Real,beamMomentTwist,0.,,"twist moment calculated for the beam.[Pa.m]"))
			((Real,beamMomentBending,0.,,"twist moment calculated for the beam.[Pa.m]"))
			((Real,initD,0.,,"equilibrium distance for interacting particles. Computed as the interparticular distance at first contact detection."))
			((bool,isBroken,false,,"flag for broken interactions"))
			((bool,isCohesive,false,,"If false, particles interact in a frictional way. If true, particles are bonded."))
			((Real,tanFrictionAngle,0.,,"tangent of Coulomb friction angle for this interaction (auto. computed). [-]"))
			((Real,beamRadius,0.,,"beamRadius=lambda*Rmin. [m]"))
			((Real,beamArea,0.,,"beamArea=pi*beamRadius^2. [m^2]"))
			((Real,beamMomInertia,0.,,"beamMomInertia=(1/4)*pi*R^4. [m^4]"))
			((Real,beamPolarMomInertia,0.,,"beamPolarMomInertia=(1/2)*pi*R^4. [m]"))
			((Real,previousDisplacement,0.,,"displacement between two particles from the previous step. [m]"))
			((bool,breakOccurred,0,,"Flag used to trigger retriangulation as soon as a cohesive bond breaks in FlowEngine (for DFNFlow use only)"))
			,
			createIndex();
			,
		);
		DECLARE_LOGGER;
		REGISTER_CLASS_INDEX(BPMPhys,NormShearPhys);
};
REGISTER_SERIALIZABLE(BPMPhys);

/** 2d functor creating InteractionPhysics (Ip2) taking BPMMat and BPMMat of 2 bodies, returning type BPMPhys */
class Ip2_BPMMat_BPMMat_BPMPhys: public IPhysFunctor{
	public:
		virtual void go(const shared_ptr<Material>& pp1, const shared_ptr<Material>& pp2, const shared_ptr<Interaction>& interaction);
		
		FUNCTOR2D(BPMMat,BPMMat);
		DECLARE_LOGGER;
		YADE_CLASS_BASE_DOC_ATTRS(Ip2_BPMMat_BPMMat_BPMPhys,IPhysFunctor,"Converts 2 :yref:`BPMMat` instances to one :yref:`BPMPhys` instance, with corresponding parameters.",                   
			((int,cohesiveTresholdIteration,1,,"should new contacts be cohesive? If strictly negativ, they will in any case. If positiv, they will before this iter, they won't afterward."))
		);
		
};
REGISTER_SERIALIZABLE(Ip2_BPMMat_BPMMat_BPMPhys);

/** 2d functor creating the interaction law (Law2) based on SphereContactGeometry (ScGeom) and BPMPhys of 2 bodies, returning type BondedContactM */
class Law2_ScGeom_BPMPhys_BondedContactM: public LawFunctor{
	public:
		virtual bool go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I);
		FUNCTOR2D(ScGeom,BPMPhys);

		YADE_CLASS_BASE_DOC_ATTRS(Law2_ScGeom_BPMPhys_BondedContactM,LawFunctor,"Interaction law for cohesive frictional material, e.g. rock that can be mechanically described with a smooth contact logic [Ivars2011].",
			((bool,neverErase,false,,"Keep interactions even if particles go away from each other (only in case another constitutive law is in the scene"))
			((bool,cracksFileExist,false,,"if true (and if :yref:`recordCracks<Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM.recordCracks>`), data are appended to an existing 'cracksKey' text file; otherwise its content is reset."))
			((string,Key,"",,"string specifying the name of saved file 'cracks___.txt', when :yref:`recordCracks<Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM.recordCracks>` is true."))
			((bool,recordCracks,false,,"if true, data about interactions that lose their cohesive feature are stored in the text file cracksKey.txt (see :yref:`Key<Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM.Key>` and :yref:`cracksFileExist<Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM.cracksFileExist>`). It contains 9 columns: the break iteration, the 3 coordinates of the contact point, the type (1 means shear break, while 0 corresponds to tensile break), the ''cross section'' (mean radius of the 2 spheres) and the 3 coordinates of the contact normal."))
			((int,nbTensCracks,0,,"number of tensile microcracks."))
			((int,nbShearCracks,0,,"number of shear microcracks."))
			((Real,totalTensCracksE,0.,,"calculate the overall energy dissipated by interparticle microcracking in tension."))
			((Real,totalShearCracksE,0.,,"calculate the overall energy dissipated by interparticle microcracking in shear."))
                        ((Real,totalCracksSurface,0.,,"calculate the total cracked surface."))
		);
		DECLARE_LOGGER;	
};
REGISTER_SERIALIZABLE(Law2_ScGeom_BPMPhys_BondedContactM);
