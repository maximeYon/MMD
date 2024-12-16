Functions for analyzing gradient waveforms

We use three basic units. First, the gradient waveform defined as a Nx3 
vector with the last three dimensions representing x, y, and z. This is
the gradient waveform played out by the gradient amplifier.

Second, we have the effect of the rf-waveform. It is a Nx1 vector with
numbers between 1 and -1. To compute the effective gradient waveform, we
need to multiply the gradient waveform played by the amplifier by the 
effect of the rf-waveform.

Third, we need to know the time step.

We refer to these variables as gwf, rf, and dt. As usual, we use SI units
for all physical entitied.
