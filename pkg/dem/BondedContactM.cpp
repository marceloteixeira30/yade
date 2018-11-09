/* MarceloDeSimone2018  */

#include"BondedContactM.hpp"
#include<core/Scene.hpp>
#include<pkg/dem/ScGeom.hpp>
#include<core/Omega.hpp>

YADE_PLUGIN((BPMMat)(BPMState)(BPMPhys)(Ip2_BPMMat_BPMMat_BPMPhys)(Law2_ScGeom_BPMPhys_BondedContactM));


/********************** Law2_ScGeom_BPMPhys_BondedContactM ****************************/
CREATE_LOGGER(Law2_ScGeom_BPMPhys_BondedContactM);

bool Law2_ScGeom_BPMPhys_BondedContactM::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact){

	const int &id1 = contact->getId1();
	const int &id2 = contact->getId2();
	ScGeom* geom = static_cast<ScGeom*>(ig.get()); 
	BPMPhys* phys = static_cast<BPMPhys*>(ip.get());
	
	Body* b1 = Body::byId(id1,scene).get();
	Body* b2 = Body::byId(id2,scene).get();

	Real Dtensile=phys->FnMax/phys->kn;
	
	string fileCracks = "cracks_"+Key+".txt";
	/// Defines the interparticular distance used for computation
	Real D = 0;

	/* This is for setting the equilibrium distance between all cohesive elements at the first contact detection*/
	if ( contact->isFresh(scene) ) { 
	  phys->normalForce = Vector3r::Zero(); 
	  phys->shearForce = Vector3r::Zero();
	  phys->initD = geom->penetrationDepth;
	}
	
	Real D = geom->penetrationDepth;
	Real cohesive_D = geom->penetrationDepth - phys->initD;

	/* Frictional contact due to overlap*/
	if (D > 0.0)
	{
	  /* NormalForce */
	  Real Fn = 0;
	  Fn = phys->kn*D;
	  
	  /* ShearForce */
	  Vector3r& shearForce = phys->shearForce;
	  shearForce = geom->rotate(phys->shearForce);
	  const Vector3r& incrementalShear = geom->shearIncrement();
	  shearForce -= phys->ks*incrementalShear;
	}
	
	if (phys->isCohesive)
	{
	  /* Normal force from the beam*/
	  Real prevD = phys->previousDisplacement
	  Real beamIncFn = beamNormalStiffness * phys->beamArea * (cohesive_D-prevD);
	  Real beamNForce = phys->beamNormalForce;
	  beamNForce = beamNForce + beamIncFn;
	  
	  /* Shear force from the beam*/
	  
	  
	  	  
	  const Real& dt = scene->dt;
	  Vector3r relAngVel = geom->getRelAngVel(de1,de2,dt);
	  /* Bending moment from the beam*/
	  Vector3r relAngVelBend = relAngVel - geom->normal.dot(relAngVel)*geom->normal; // keep only the bending part
	  Vector3r relRotBend = relAngVelBend*dt; // relative rotation due to rolling behaviour	
	  // incremental formulation for the bending moment (as for the shear part)
	  Vector3r& momentBend = phys->beamMomentBending;
	  momentBend = geom->rotate(momentBend); // rotate moment vector (updated)
	  momentBend = momentBend-phys->kr*relRotBend;
	  
	  // *** Torsion ***//
	  Vector3r relAngVelTwist = geom->normal.dot(relAngVel)*geom->normal;
	  Vector3r relRotTwist = relAngVelTwist*dt; // component of relative rotation along n  FIXME: sign?
	  // incremental formulation for the torsional moment
	  Vector3r& momentTwist = phys->beamMomentTwist;
	  momentTwist = geom->rotate(momentTwist); // rotate moment vector (updated)
	  momentTwist = momentTwist-phys->ktw*relRotTwist; // FIXME: sign?
	  
	}
	
	/* Apply forces */
	phys->normalForce = Fn*geom->normal;
	
	Vector3r f = phys->normalForce + shearForce;
	
	/// applyForceAtContactPoint computes torque also and, for now, we don't want rotation for particles on joint (some errors in calculation due to specific geometry) 
 	//applyForceAtContactPoint(f, geom->contactPoint, I->getId2(), b2->state->pos, I->getId1(), b1->state->pos, scene);
	scene->forces.addForce (id1,-f);
	scene->forces.addForce (id2, f);
	
	// simple solution to avoid torque computation for particles interacting on a smooth joint 
	if ( (phys->isOnJoint)&&(smoothJoint) ) return true;
	
	/// those lines are needed if rootBody->forces.addForce and rootBody->forces.addMoment are used instead of applyForceAtContactPoint -> NOTE need to check for accuracy!!!
	scene->forces.addTorque(id1,(geom->radius1-0.5*geom->penetrationDepth)* geom->normal.cross(-f));
	scene->forces.addTorque(id2,(geom->radius2-0.5*geom->penetrationDepth)* geom->normal.cross(-f));
	return true;
	
}

CREATE_LOGGER(Ip2_BPMMat_BPMMat_BPMPhys);

void Ip2_BPMMat_BPMMat_BPMPhys::go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction){

	/* avoid updates of interaction if it already exists */
	if( interaction->phys ) return; 

	ScGeom* geom=dynamic_cast<ScGeom*>(interaction->geom.get());
	assert(geom);

	const shared_ptr<BPMMat>& yade1 = YADE_PTR_CAST<BPMMat>(b1);
	const shared_ptr<BPMMat>& yade2 = YADE_PTR_CAST<BPMMat>(b2);
	BPMState* st1=dynamic_cast<BPMState*>(Body::byId(interaction->getId1(),scene)->state.get());
	BPMState* st2=dynamic_cast<BPMState*>(Body::byId(interaction->getId2(),scene)->state.get());
	
	shared_ptr<BPMPhys> contactPhysics(new BPMPhys());
	
	/* From material properties */
	Real E1 	= yade1->young;
	Real E2 	= yade2->young;
	Real v1 	= yade1->poisson;
	Real v2 	= yade2->poisson;
	Real f1 	= yade1->frictionAngle;
	Real f2 	= yade2->frictionAngle;
	Real rf1 	= yade1->residualFrictionAngle>=0? yade1->residualFrictionAngle: yade1->frictionAngle;
	Real rf2 	= yade2->residualFrictionAngle>=0? yade2->residualFrictionAngle: yade2->frictionAngle;
	Real norCoh1	= yade1->normalCohesion;
	Real norCoh2	= yade2->normalCohesion;
	Real sheCoh1	= yade1->shearCohesion;
	Real sheCoh2	= yade2->shearCohesion;
	Real lambda1	= yade1->lambda;
	Real lambda2	= yade2->lambda;

	/* From interaction geometry */
	Real R1= geom->radius1;
	Real R2= geom->radius2;
	
	/* From contact properties*/
	contactPhysics->beamRadius = min(lambda1,lambda2)*pow(min(R1,R2),2);
	contactPhysics->beamArea = Mathr::PI*contactPhysics->beamRadius;
	contactPhysics->beamMomInertia = (1/4)*Mathr::PI*pow(contactPhysics->beamRadius,4);
	contactPhysics->beamPolarMomInertia = (1/2)*Mathr::PI*pow(contactPhysics->beamRadius,4);

	/* Pass values to BPMPhys.*/
	
	// elastic properties (ks and kn)
	contactPhysics->kn = 2.*E1*R1*E2*R2/(E1*R1+E2*R2);
        ( (v1==0)&&(v2==0) )? contactPhysics->ks=0 : contactPhysics->ks = 2.*E1*R1*v1*E2*R2*v2/(E1*R1*v1+E2*R2*v2);
	
	// cohesive properties
	///to set if the contact is cohesive or not
	if ( ((cohesiveTresholdIteration < 0) || (scene->iter < cohesiveTresholdIteration)) && (std::min(SigT1,SigT2)>0 || std::min(Coh1,Coh2)>0) && (yade1->type == yade2->type)){ 
	  contactPhysics->isCohesive=true;
	  st1->nbInitBonds++;
	  st2->nbInitBonds++;
	}
	
	if ( contactPhysics->isCohesive ) {
	  contactPhysics->FnMax = std::min(SigT1,SigT2)*contactPhysics->crossSection;
	  contactPhysics->FsMax = std::min(Coh1,Coh2)*contactPhysics->crossSection;
	}
	
        // frictional properties      
        contactPhysics->isCohesive? contactPhysics->tanFrictionAngle = std::tan(std::min(f1,f2)) : contactPhysics->tanFrictionAngle = std::tan(std::min(rf1,rf2));

	interaction->phys = contactPhysics;
}

BPMPhys::~BPMPhys(){}
