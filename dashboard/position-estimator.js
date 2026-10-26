(function(root){
"use strict";
const G=9.81,DEG=Math.PI/180;
function norm(v){return Math.hypot(v[0],v[1],v[2])}
function rotation(e){
  const r=e[0]*DEG,p=e[1]*DEG,y=e[2]*DEG,cr=Math.cos(r),sr=Math.sin(r),cp=Math.cos(p),sp=Math.sin(p),cy=Math.cos(y),sy=Math.sin(y);
  return[cy*cp,cy*sp*sr-sy*cr,cy*sp*cr+sy*sr,sy*cp,sy*sp*sr+cy*cr,sy*sp*cr-cy*sr,-sp,cp*sr,cp*cr];
}
function mul(m,v){return[m[0]*v[0]+m[1]*v[1]+m[2]*v[2],m[3]*v[0]+m[4]*v[1]+m[5]*v[2],m[6]*v[0]+m[7]*v[1]+m[8]*v[2]]}
class Butterworth{
  constructor(cutoffHz,highPass){this.cutoffHz=cutoffHz;this.highPass=highPass;this.reset()}
  reset(){this.initialized=false;this.state=0}
  process(value,dt){
    if(!(dt>0))return 0;
    const normalized=Math.min(this.cutoffHz*dt,.49),k=Math.tan(Math.PI*normalized),divisor=1+k,b0=(this.highPass?1:k)/divisor,b1=(this.highPass?-1:k)/divisor,a1=(k-1)/divisor;
    if(!this.initialized){this.state=(b1-a1*b0)/(1+a1)*value;this.initialized=true}
    const output=b0*value+this.state;this.state=b1*value-a1*output;return output;
  }
}
class PositionEstimator{
  constructor(options={}){
    this.config=Object.assign({nominalRate:800,maximumDeltaSeconds:.05,stationaryThresholdMps2:.08,stationaryHighPassHz:.001,stationaryLowPassHz:5,stationaryConfirmationSeconds:.075,motionPreRollSeconds:.04,initialCalibrationSeconds:.35,staticCalibrationTauSeconds:10,historySeconds:25},options);
    /* Backward-compatible test/UI option name. */
    if(options.calibrationSeconds!==undefined)this.config.initialCalibrationSeconds=options.calibrationSeconds;
    this.reset();
  }
  reset(){
    this.position=[0,0,0];this.velocity=[0,0,0];this.linearAcceleration=[0,0,0];this.lastSeq=null;this.rate=this.config.nominalRate;this.rateAnchor=null;this.clock=0;this.history=[];this.segmentStart=null;this.restCandidateStart=null;this.restCandidateSeconds=0;this.isStatic=true;this.distance=0;this.revision=0;this.uncertainty=.001;this.filteredMotion=0;this.accelerationReferenceMps2=G;this.calibrationSum=0;this.calibrationSamples=0;this.calibrationSeconds=0;this.calibrated=false;this.highPass=new Butterworth(this.config.stationaryHighPassHz,true);this.lowPass=new Butterworth(this.config.stationaryLowPassHz,false);
  }
  resetOrigin(){this.position=[0,0,0];this.velocity=[0,0,0];this.linearAcceleration=[0,0,0];this.history=[];this.segmentStart=null;this.restCandidateStart=null;this.restCandidateSeconds=0;this.isStatic=true;this.distance=0;this.clock=0;this.uncertainty=.001;this.revision++}
  updateRate(seq,hostMs){
    if(!this.rateAnchor){this.rateAnchor={seq,hostMs};return}
    const ds=(seq-this.rateAnchor.seq)&65535,elapsed=(hostMs-this.rateAnchor.hostMs)/1000;
    if(elapsed>=.4&&ds>0){const measured=ds/elapsed,n=this.config.nominalRate;if(measured>n*.75&&measured<n*1.25)this.rate=this.rate*.85+measured*.15;this.rateAnchor={seq,hostMs}}
  }
  process(sample){
    const seq=sample.seq&65535,hostMs=sample.hostMs,requestedRate=sample.nominalRate||this.config.nominalRate;
    if(requestedRate!==this.config.nominalRate){this.config.nominalRate=requestedRate;this.rate=requestedRate;this.rateAnchor=null}
    this.updateRate(seq,hostMs);
    let ds=this.lastSeq===null?0:(seq-this.lastSeq)&65535;if(ds>32)ds=0;
    let dt=ds>0?ds/Math.max(this.rate,1):0;if(dt>this.config.maximumDeltaSeconds)dt=0;this.lastSeq=seq;
    const rawBody=sample.accG.map(v=>v*G),rawMagnitude=norm(rawBody),high=Math.abs(this.highPass.process(rawMagnitude,dt));
    this.filteredMotion=this.lowPass.process(high,dt);
    const below=dt===0||this.filteredMotion<this.config.stationaryThresholdMps2,wasStatic=this.isStatic;
    let staticNow=wasStatic?below:false,completedSegmentEnd=null;
    if(!wasStatic){
      if(below){
        if(this.restCandidateStart===null)this.restCandidateStart=this.history.length;
        this.restCandidateSeconds+=dt;
        if(this.restCandidateSeconds>=this.config.stationaryConfirmationSeconds){staticNow=true;completedSegmentEnd=this.restCandidateStart}
      }else{this.restCandidateStart=null;this.restCandidateSeconds=0}
    }
    if(wasStatic&&!staticNow){
      const pre=Math.max(0,Math.round(this.config.motionPreRollSeconds/Math.max(dt,1e-6)));
      this.segmentStart=Math.max(0,this.history.length-pre);this.replayMotionStart(this.segmentStart);this.restCandidateStart=null;this.restCandidateSeconds=0;
    }
    if(staticNow&&rawMagnitude>1e-9){
      if(!this.calibrated){
        this.calibrationSum+=rawMagnitude;this.calibrationSamples++;this.calibrationSeconds+=dt;this.accelerationReferenceMps2=this.calibrationSum/this.calibrationSamples;this.calibrated=this.calibrationSeconds>=this.config.initialCalibrationSeconds;
      }else if(dt>0){const a=1-Math.exp(-dt/this.config.staticCalibrationTauSeconds);this.accelerationReferenceMps2+=a*(rawMagnitude-this.accelerationReferenceMps2)}
    }
    const scale=G/Math.max(this.accelerationReferenceMps2,1e-9),world=mul(rotation(sample.eulerDeg),rawBody.map(v=>v*scale));
    this.linearAcceleration=[world[0],world[1],world[2]-G];
    if(!this.calibrated||dt===0)this.linearAcceleration=[0,0,0];
    const previousP=this.position.slice(),previousV=this.velocity.slice(),rawV=previousV.map((v,i)=>v+this.linearAcceleration[i]*dt),nextV=staticNow?[0,0,0]:rawV.slice(),nextP=staticNow?previousP:previousP.map((p,i)=>p+nextV[i]*dt);
    if(!staticNow)this.distance+=norm(nextV)*dt;
    this.clock+=dt;this.history.push({time:this.clock,position:nextP,velocity:nextV,rawVelocity:rawV,linearAcceleration:this.linearAcceleration.slice(),dt,isStatic:staticNow});this.position=nextP;this.velocity=nextV;this.isStatic=staticNow;
    if(!wasStatic&&staticNow&&this.segmentStart!==null){this.correctSegment(this.segmentStart,completedSegmentEnd);this.segmentStart=null;this.restCandidateStart=null;this.restCandidateSeconds=0}
    this.trim();this.uncertainty=staticNow?Math.max(.001,this.uncertainty*.995):Math.min(9.999,this.uncertainty+dt*(.008+.06*norm(this.linearAcceleration)));
    return this.output();
  }
  replayMotionStart(start){
    let p=start?this.history[start-1].position.slice():[0,0,0],v=[0,0,0];
    for(let i=start;i<this.history.length;i++){const r=this.history[i];r.rawVelocity=v.map((x,a)=>x+r.linearAcceleration[a]*r.dt);r.velocity=r.rawVelocity.slice();r.position=p.map((x,a)=>x+r.velocity[a]*r.dt);r.isStatic=false;p=r.position.slice();v=r.velocity.slice()}
    this.position=p;this.velocity=v;this.revision++;
  }
  correctSegment(start,end){
    const count=end-start;if(count<=0)return;const endV=this.history[end-1].rawVelocity.slice();let p=start?this.history[start-1].position.slice():[0,0,0];
    for(let i=start;i<end;i++){const r=this.history[i],q=(i-start+1)/count;r.velocity=r.rawVelocity.map((v,a)=>v-endV[a]*q);r.position=p.map((x,a)=>x+r.velocity[a]*r.dt);p=r.position.slice()}
    for(let i=end;i<this.history.length;i++){this.history[i].position=p.slice();this.history[i].velocity=[0,0,0];this.history[i].isStatic=true}
    this.position=p;this.velocity=[0,0,0];this.revision++;
  }
  trim(){
    if(this.segmentStart!==null)return;const cutoff=this.clock-this.config.historySeconds;let n=0;
    while(n+2<this.history.length&&this.history[n+1].time<cutoff)n++;
    if(n)this.history.splice(0,n);
  }
  points(seconds=20){const cutoff=this.clock-seconds,last=[];let lastTime=-Infinity;for(const r of this.history)if(r.time>=cutoff&&r.time-lastTime>=.025){last.push({time:r.time,position:r.position.slice(),isStatic:r.isStatic});lastTime=r.time}return last}
  output(){return{position:this.position.slice(),velocity:this.velocity.slice(),linearAcceleration:this.linearAcceleration.slice(),distance:this.distance,uncertainty:this.uncertainty,isStatic:this.isStatic,deviceStatic:null,calibrated:this.calibrated,calibrationProgress:Math.min(1,this.calibrationSeconds/this.config.initialCalibrationSeconds),rate:this.rate,revision:this.revision,filteredMotion:this.filteredMotion,accelerationReferenceMps2:this.accelerationReferenceMps2}}
}
root.CJ02Position={PositionEstimator,rotation,mul};
})(typeof window!=="undefined"?window:globalThis);
