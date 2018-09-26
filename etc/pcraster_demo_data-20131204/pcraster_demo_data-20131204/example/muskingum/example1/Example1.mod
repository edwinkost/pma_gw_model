# example from Applied Hydrology book
timer
  1 40 1;
initial
  TimeStepInSecs = 3600;
  K = spatial(scalar(2.3))/TimeStepInSecs;  # [sec]
  x = spatial(0.15);
  PrevFlowRate = if(id.map eq 1 then 85 else 85)/TimeStepInSecs ;  # m3/sec
  PrevInflow = if(id.map eq 1 then 85 else 0)/TimeStepInSecs;      # m3/sec
  SegmentLength = scalar(1);
  Iterations = nominal(1);
dynamic
  InflowTimestep = timeinputscalar(inflowTimestep.tss, id.map)/TimeStepInSecs;  # m3/sec
  FlowRate = muskingum(ldd.map,PrevFlowRate,PrevInflow,InflowTimestep,K,x,SegmentLength,Iterations,TimeStepInSecs);
report flowRateTimestep.tss = timeoutput(id.map,FlowRate*TimeStepInSecs);
  PrevFlowRate = FlowRate;
  PrevInflow = InflowTimestep;
