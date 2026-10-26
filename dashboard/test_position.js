const fs=require("fs"),vm=require("vm"),assert=require("assert"),ctx={console,Math};
ctx.globalThis=ctx;vm.createContext(ctx);
vm.runInContext(fs.readFileSync(__dirname+"/position-estimator.js","utf8"),ctx);
const {PositionEstimator:E,rotation:R}=ctx.CJ02Position,G=9.80665;
const near=(actual,expected,tolerance,message)=>assert(Math.abs(actual-expected)<=tolerance,`${message}: ${actual} vs ${expected}`);
const transpose=(m,v)=>[m[0]*v[0]+m[3]*v[1]+m[6]*v[2],m[1]*v[0]+m[4]*v[1]+m[7]*v[2],m[2]*v[0]+m[5]*v[1]+m[8]*v[2]];

function rig(rate=100,euler=[0,0,0],bias=[0,0,0]){
  const est=new E({nominalRate:rate,calibrationSeconds:.5,historySeconds:25});let seq=0,hostMs=0,out;
  function feed(seconds,worldLinear=[0,0,0],isStatic=false){
    const samples=Math.round(seconds*rate),matrix=R(euler);
    for(let i=0;i<samples;i++){
      const body=transpose(matrix,[worldLinear[0],worldLinear[1],worldLinear[2]+G]).map((v,a)=>(v+bias[a])/G);
      out=est.process({seq:seq++&65535,hostMs,accG:body,eulerDeg:euler,isStatic,nominalRate:rate});hostMs+=1000/rate;
    }
    return out;
  }
  return{est,feed,get out(){return out}};
}

for(const euler of [[0,0,0],[30,-20,45],[-70,35,-120]]){
  const q=rig(200,euler);q.feed(.2,[0,0,0],true);
  assert(!q.out.calibrated,"calibration completed too early");
  q.feed(.5,[0,0,0],true);assert(q.out.calibrated);
  q.out.linearAcceleration.forEach((v,i)=>near(v,0,2e-3,`tilted gravity axis ${i}`));
}

const motion=rig(100,[25,-15,40]);motion.feed(.7,[0,0,0],true);
motion.feed(1,[2,0,0],false);motion.feed(1,[-2,0,0],false);motion.feed(.2,[0,0,0],true);
near(motion.out.velocity[0],0,1e-9,"ZUPT velocity");near(motion.out.position[0],2,.08,"two metre trajectory");
near(motion.out.position[1],0,2e-3,"cross-axis Y");near(motion.out.position[2],0,2e-3,"cross-axis Z");

const longMove=rig(100);longMove.feed(.7,[0,0,0],true);longMove.feed(15,[2,0,0],false);longMove.feed(15,[-2,0,0],false);longMove.feed(.2,[0,0,0],true);
assert(longMove.est.history.length>2500,"active segment was incorrectly trimmed");near(longMove.out.position[0],450,2,"long segment position");

const switched=rig(800);switched.feed(.7,[0,0,0],true);let switchSeq=1000,switchHost=1000;
for(const rate of [1600,800,400,200]){
  const out=switched.est.process({seq:switchSeq++&65535,hostMs:switchHost,accG:[0,0,1],eulerDeg:[0,0,0],isStatic:true,nominalRate:rate});
  near(out.rate,rate,1e-9,`rate switch ${rate}`);switchHost+=1000/rate;
}

motion.est.resetOrigin();assert.deepStrictEqual(Array.from(motion.est.output().position),[0,0,0]);
console.log("PASS: reference calibration, tilted gravity, 3D integration, ZUPT, long segment, 4 rates, reset");

if(process.argv[2]){
  const path=require("path"),dir=process.argv[2],raw=fs.readFileSync(path.join(dir,"raw.csv"),"utf8").trim().split(/\r?\n/).slice(1),att=fs.readFileSync(path.join(dir,"attitude.csv"),"utf8").trim().split(/\r?\n/).slice(1),n=Math.min(raw.length,att.length),est=new E({nominalRate:800,calibrationSeconds:.5});let out,maxStaticLinear=0,accNormSum=0;
  for(let i=0;i<n;i++){const r=raw[i].split(",").map(Number),a=att[i].split(",").map(Number),accG=r.slice(2,5).map(x=>x/2048);accNormSum+=Math.hypot(...accG);out=est.process({seq:r[1],hostMs:r[0]*1000,accG,eulerDeg:a.slice(2,5),isStatic:Boolean(a[6]&4),nominalRate:800});if(out.calibrated&&out.isStatic)maxStaticLinear=Math.max(maxStaticLinear,Math.hypot(...out.linearAcceleration))}
  const meanG=accNormSum/n;assert(meanG>.5&&meanG<1.5,`recording is not current +/-16 g data: mean=${meanG.toFixed(3)} g`);
  assert(out.calibrated);assert(out.position.every(Number.isFinite));assert(out.velocity.every(Number.isFinite));
  console.log(`PASS: recorded hardware ${n} frames, final=${out.position.map(x=>x.toFixed(4)).join(",")} m, max-static-linear=${maxStaticLinear.toFixed(4)} m/s2, rate=${out.rate.toFixed(1)} Hz`);
}
