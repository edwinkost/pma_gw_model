binding
  DEM = topo.map;
timer
  1 20 1;
initial
  TimeStepInSecs = 3600;
  K = spatial(scalar(2.3))/TimeStepInSecs;  # [sec]
  x = spatial(0.15);
  SegmentLength = spatial(scalar(1));
  Iterations = nominal(10);
  lddMap = lddcreate(DEM,1e37,1e37,1e37,1e37);
  ups = accuflux(lddMap,1);
  pitNumber = if(ups eq mapmaximum(ups) then nominal(1) else 0);
  Channel = ups gt 100;
  LocalBaseFlow = 0.01;
  lddChannel = if(Channel then lddMap);
  lddOverlandFlow = if(Channel then 5 else lddMap);

  PrevFlowRate = iniFlowRate.map;
  PrevInflow = accuflux(lddOverlandFlow,LocalBaseFlow)/TimeStepInSecs;

dynamic
  LocalSurfWater = timeinputscalar(inflowTimestep.tss, 1);  # m3/timestep
  overlandFlow = accuflux(lddOverlandFlow,LocalBaseFlow+LocalSurfWater); #m3/timestep

  inFlow= overlandFlow/TimeStepInSecs;  #m3/sec 

  report FlowRate = muskingum(lddMap,PrevFlowRate,PrevInflow,
                        inFlow,K,x,
                        SegmentLength, Iterations,TimeStepInSecs);

  report flowRate.tss = timeoutput(pitNumber,FlowRate*TimeStepInSecs);
  PrevFlowRate = FlowRate;  # m3/sec
  PrevInflow = inFlow;      # m3/sec
