/* LucScholtes2010  */

#include"BPM.hpp"
#include<core/Scene.hpp>
#include<pkg/dem/ScGeom.hpp>
#include<core/Omega.hpp>

YADE_PLUGIN((BPMpmMat)(BPMpmState)(BPMpmPhys)(Ip2_BPMpmMat_BPMpmMat_BPMpmPhys)(Law2_ScGeom_BPMpmPhys_BondedContactM));


/********************** Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM ****************************/
CREATE_LOGGER(Law2_ScGeom_BPMpmPhys_BondedContactM);

bool Law2_ScGeom_BPMpmPhys_BondedContactM::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact){

	const int &id1 = contact->getId1();
	const int &id2 = contact->getId2();
	ScGeom* geom = static_cast<ScGeom*>(ig.get()); 
	BPMpmPhys* phys = static_cast<BPMpmPhys*>(ip.get());
	
	Body* b1 = Body::byId(id1,scene).get();
	Body* b2 = Body::byId(id2,scene).get();
	
	string fileCracks = "cracks_"+Key+".txt";
	/// Defines the interparticular distance used for computation
	Real D = 0;

	/*this is for setting the equilibrium distance between all cohesive elements at the first contact detection*/
	if ( contact->isFresh(scene) ) { 
	  phys->normalForce = Vector3r::Zero(); 
	  phys->shearForce = Vector3r::Zero();
	  phys->initD = geom->penetrationDepth; 
	}
	
	D = geom->penetrationDepth;
	Real cohesive_D = geom->penetrationDepth - phys->initD;
	
	/* Determination of interaction */
	if (D < 0) { //tensile configuration
	  if ( !phys->isCohesive) {
	    if (!neverErase) return false;
	    else {
	      phys->shearForce = Vector3r::Zero();
	      phys->normalForce = Vector3r::Zero();
	      phys->isCohesive =0;
	      return true;
	    }
	  }
	}
	
	Vector3r& shearForce = phys->shearForce;
	const Vector3r& incrementalShear = geom->shearIncrement();
	
	// Cohesive variables declaration
	Vector3r& beamNForce = phys->beamNormalForce;
	Vector3r& beamSForce = phys->beamShearForce;
	Vector3r& momentBend = phys->beamMomentBending;
	Vector3r& momentTwist = phys->beamMomentTwist;
	Real normalStress = 0.;
	Real shearStress = 0.;
	
	if ( phys->isCohesive ) {
	  /* Normal force from the beam */
	  Real beamIncFn = phys->beamNormalStiffness * phys->beamArea * (cohesive_D-phys->previousDisplacement);
	  Vector3r beamIncForce = beamIncFn * geom->normal;
	  beamNForce = beamNForce + beamIncForce;
	  
	  /* Shear force from the beam */
	  beamSForce = geom->rotate(phys->beamShearForce);
	  beamSForce -= phys->beamShearStiffness * phys->beamArea * incrementalShear;
	  
	  /* Moments from the beam */
	  const Real& dt = scene->dt;
	  State* de1 = Body::byId(id1,scene)->state.get();
	  State* de2 = Body::byId(id2,scene)->state.get();
	  Vector3r relAngVel = geom->getRelAngVel(de1,de2,dt);
	  
	  /* Bending moment */
	  Vector3r relAngVelBend = relAngVel - geom->normal.dot(relAngVel)*geom->normal; // keep only the bending part
	  Vector3r relRotBend = relAngVelBend*dt; // relative rotation due to rolling behaviour	
	  // incremental formulation for the bending moment (as for the shear part)
	  momentBend = geom->rotate(momentBend);
	  momentBend -= phys->beamNormalStiffness*phys->beamMomInertia*relRotBend;
	  
	  /* Torsion */
	  Vector3r relAngVelTwist = geom->normal.dot(relAngVel)*geom->normal;
	  Vector3r relRotTwist = relAngVelTwist*dt;
	  // incremental formulation for the torsional moment
	  momentTwist = geom->rotate(momentTwist);
	  momentTwist -= phys->beamShearStiffness*phys->beamPolarMomInertia*relRotTwist;
	  
	  /* Limits for the moments and forces */
	  Real sign = (cohesive_D > 0.0) ? 1.0 : ((cohesive_D < 0.0) ? -1.0 : 0.0);
	  normalStress = -((sign * beamNForce.norm()) / phys->beamArea) + (phys->beamBeta * (momentBend.norm() * phys->beamRadius) / phys->beamMomInertia);
	  shearStress = (beamSForce.norm() / phys->beamArea) + (phys->beamBeta * (momentTwist.norm() * phys->beamRadius) / phys->beamPolarMomInertia);
	  
	  /* Update shear resistance */
	  Real& beamSCoh = phys->beamSCoh;
	  beamSCoh = (phys->beamShearCohesion + (sign * beamNForce.norm() / phys->beamArea) * phys->tanFrictionAngle);
	  
	  if (fabs(shearStress) > fabs(beamSCoh))
	  {
	    nbShearCracks++;
	    phys->isCohesive = 0;
	    phys->breakOccurred = true;  // flag to trigger remesh for DFNFlowEngine
	    phys->isBroken = true; // flag for DFNFlowEngine
	    
	    // update body state with the number of broken bonds
	    BPMpmState* st1=dynamic_cast<BPMpmState*>(b1->state.get());
	    BPMpmState* st2=dynamic_cast<BPMpmState*>(b2->state.get());
	    st1->nbBrokenBonds++;
	    st2->nbBrokenBonds++;
	    st1->damageIndex+=1.0/st1->nbInitBonds;
	    st2->damageIndex+=1.0/st2->nbInitBonds;
	    
	    if ( D < 0 ) { // spheres do not touch
                if (!neverErase) return false;
                else {
                    phys->shearForce = Vector3r::Zero();
                    phys->normalForce = Vector3r::Zero();
                    return true;
                }
	    }
	  }
	  
	  if (fabs(normalStress) > fabs(phys->beamNormalCohesion) && cohesive_D < 0.0)
	  {
	    nbTensCracks++;
	    phys->isCohesive = 0;
	    phys->breakOccurred = true;  // flag to trigger remesh for DFNFlowEngine
	    phys->isBroken = true; // flag for DFNFlowEngine
	    
            // update body state with the number of broken bonds
	    BPMpmState* st1=dynamic_cast<BPMpmState*>(b1->state.get());
	    BPMpmState* st2=dynamic_cast<BPMpmState*>(b2->state.get());
            st1->nbBrokenBonds++;
	    st2->nbBrokenBonds++;
	    st1->damageIndex+=1.0/st1->nbInitBonds;
	    st2->damageIndex+=1.0/st2->nbInitBonds;
	    
	    if (!neverErase) return false; 
	    else {
	      phys->shearForce = Vector3r::Zero();
	      phys->normalForce = Vector3r::Zero();
	      return true;
	    }
	  }
	}
	
	Real Fn = 0;
	if (D > 0 && cohesive_D > 0)
	{
	  /* NormalForce */
	  Fn = phys->kn*cohesive_D; 
        
	  /* ShearForce */
	  shearForce = geom->rotate(phys->shearForce);
	  shearForce -= phys->ks*incrementalShear;
	
	  /* Mohr-Coulomb criterion */
	  Real maxFs = Fn*phys->tanResidualFrictionAngle;
	  Real scalarShearForce = shearForce.norm();
               
	  if (scalarShearForce > maxFs) {
	    if (scalarShearForce != 0)
	      shearForce*=maxFs/scalarShearForce;
	    else
	      shearForce=Vector3r::Zero();
	  }
	}
	
	/* Apply forces */
	phys->normalForce = Fn*geom->normal;
	
	Vector3r f = phys->normalForce + beamNForce + shearForce + beamSForce;
	
	/// applyForceAtContactPoint computes torque also and, for now, we don't want rotation for particles on joint (some errors in calculation due to specific geometry) 
 	//applyForceAtContactPoint(f, geom->contactPoint, I->getId2(), b2->state->pos, I->getId1(), b1->state->pos, scene);
	scene->forces.addForce (id1,-f);
	scene->forces.addForce (id2, f);
	
	/// those lines are needed if rootBody->forces.addForce and rootBody->forces.addMoment are used instead of applyForceAtContactPoint -> NOTE need to check for accuracy!!!
	scene->forces.addTorque(id1,(geom->radius1-0.5*geom->penetrationDepth)* geom->normal.cross(-f));
	scene->forces.addTorque(id2,(geom->radius2-0.5*geom->penetrationDepth)* geom->normal.cross(-f));
	
	Vector3r beamMoment = momentBend + momentTwist;
	scene->forces.addTorque(id1,- beamMoment);
	scene->forces.addTorque(id2,beamMoment);
	
	phys->previousDisplacement = cohesive_D;

	return true;
	
}

CREATE_LOGGER(Ip2_BPMpmMat_BPMpmMat_BPMpmPhys);

void Ip2_BPMpmMat_BPMpmMat_BPMpmPhys::go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction){

	/* avoid updates of interaction if it already exists */
	if( interaction->phys ) return; 

	ScGeom* geom=dynamic_cast<ScGeom*>(interaction->geom.get());
	assert(geom);

	const shared_ptr<BPMpmMat>& yade1 = YADE_PTR_CAST<BPMpmMat>(b1);
	const shared_ptr<BPMpmMat>& yade2 = YADE_PTR_CAST<BPMpmMat>(b2);
	BPMpmState* st1=dynamic_cast<BPMpmState*>(Body::byId(interaction->getId1(),scene)->state.get());
	BPMpmState* st2=dynamic_cast<BPMpmState*>(Body::byId(interaction->getId2(),scene)->state.get());
	
	shared_ptr<BPMpmPhys> contactPhysics(new BPMpmPhys()); 
	
	/* From material properties */
	Real E1 	= yade1->young;
	Real E2 	= yade2->young;
	Real v1 	= yade1->poisson;
	Real v2 	= yade2->poisson;
	Real f1 	= yade1->frictionAngle;
	Real f2 	= yade2->frictionAngle;
        Real rf1 	= yade1->residualFrictionAngle>=0? yade1->residualFrictionAngle: yade1->frictionAngle;
	Real rf2 	= yade2->residualFrictionAngle>=0? yade2->residualFrictionAngle: yade2->frictionAngle;
	Real cohesiveE1	= yade1->cohesiveYoung;
	Real cohesiveE2	= yade2->cohesiveYoung;
	Real cohesiveV1 = yade1->cohesivePoisson;
	Real cohesiveV2 = yade2->cohesivePoisson;
	Real norCoh1	= yade1->normalCohesion;
	Real norCoh2	= yade2->normalCohesion;
	Real sheCoh1	= yade1->shearCohesion;
	Real sheCoh2	= yade2->shearCohesion;
	Real lambda1	= yade1->lambdaCohesion;
	Real lambda2	= yade2->lambdaCohesion;
	Real beta1	= yade1->betaCohesion;
	Real beta2	= yade2->betaCohesion;

	/* From interaction geometry */
	Real R1= geom->radius1;
	Real R2= geom->radius2;

	/* Pass values to JCFpmPhys. In case of a "jointed" interaction, the following values will be replaced by other ones later (in few if(){} blocks)*/
	
	// elastic properties
	contactPhysics->kn = 2.*E1*R1*E2*R2/(E1*R1+E2*R2);
        ( (v1==0)&&(v2==0) )? contactPhysics->ks=0 : contactPhysics->ks = 2.*E1*R1*v1*E2*R2*v2/(E1*R1*v1+E2*R2*v2);
	
	// cohesive properties
	///to set if the contact is cohesive or not
	if ( ((cohesiveTresholdIteration < 0) || (scene->iter < cohesiveTresholdIteration)) && (std::min(norCoh1,norCoh2)>0 || std::min(sheCoh1,sheCoh2)>0) && (yade1->type == yade2->type)){ 
	  contactPhysics->isCohesive=true;
	  st1->nbInitBonds++;
	  st2->nbInitBonds++;
	}
	
	if ( contactPhysics->isCohesive ) {
	  contactPhysics->beamRadius = min(lambda1,lambda2)*pow(min(R1,R2),2);
	  contactPhysics->beamArea = Mathr::PI*contactPhysics->beamRadius;
	  contactPhysics->beamMomInertia = (1./4.)*Mathr::PI*pow(contactPhysics->beamRadius,4);
	  contactPhysics->beamPolarMomInertia = (1./2.)*Mathr::PI*pow(contactPhysics->beamRadius,4);
	  contactPhysics->beamNormalStiffness = 2.*cohesiveE1*R1*cohesiveE2*R2/(cohesiveE1*R1+cohesiveE2*R2);
	  contactPhysics->beamShearStiffness = 2.*cohesiveE1*R1*cohesiveV1*cohesiveE2*R2*cohesiveV2/(cohesiveE1*R1*cohesiveV1+cohesiveE2*R2*cohesiveV2);
	  contactPhysics->beamNormalCohesion = std::min(norCoh1,norCoh2);
	  contactPhysics->beamShearCohesion = std::min(sheCoh1,sheCoh2);
	  contactPhysics->beamSCoh = contactPhysics->beamShearCohesion;
	  contactPhysics->beamBeta = std::min(beta1,beta2);
	}

        // frictional properties      
        contactPhysics->tanFrictionAngle = std::tan(std::min(f1,f2));
	contactPhysics->tanResidualFrictionAngle = std::tan(std::min(rf1,rf2));

	interaction->phys = contactPhysics;
}

BPMpmPhys::~BPMpmPhys(){}
