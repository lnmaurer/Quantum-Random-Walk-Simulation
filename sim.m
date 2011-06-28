%The big idea:
%To keep the code modular, resonators and QBs are handled seperately for as long as possible.
%For example, there are seperate creation/annihilation operators for QBs and resonators.
%Later, these get combined with the Kronecker tensor product to form operators for the whole system.

%The first section of the code is a bunch of functions that do useful things (like evolving states).
%The second section uses that code to simulate experiments.
%The third section is for temporary functions.

%NOTES:
%-states are column vectors
%-expected values (EVs) and straight up fock states are row vectors
%-when I say 'state' I mean state of the full system, 'resonator state' and 'QB state' refer to the states of those parts of the system
%-density matrices are stored in cell arrays
%-for operators that act on states of the full system, the correct order is kron(resonator operator, QB operator)

%constants
global N_ph = 2; %consider up to this many photons
global N_res = 2; %number of resonators to consider
global h_bar = 1.05457148e-34;
global w0 = 2*pi*5e9; %resonator resonance frequency 
global wQB = 2*pi*5e9; %QB resonance frequency
global wInt = 2*pi*1/10e-9/4; %interaction between QB and resonators
global Wcoup = 2*pi*1/100e-9/4; %coupling between resonators
global resRelaxationTime = 600e-9; %resonator relaxtion rate A RELAXATION RATE OF <=0 TURNS OF RELAXATION
global QBRelaxationTime = 200e-9; %resonator relaxtion rate A RELAXATION RATE OF <=0 TURNS OF RELAXATION
global N_resonator_states; %run generate_resonator_states to populate
global resonator_states; %run generate_resonator_states to populate
global RelEvlStepSz = 10e-12; %step size for evolution with relaxation NOT USED BY A FUNCTIONS YET

%%%%%%%%%%%%%%%%%%% USEFUL FUNCTIONS %%%%%%%%%%%%%%%%%%%

function ev = state_expectation_value(state, operator)
  ev = abs(state'*(operator*state)); %use 'abs' because small small imag part often remains due to rounding errors
end

function ev = DM_expectation_value(DM, operator)
  ev = abs(trace(DM*operator)); %use 'abs' because small small imag part often remains due to rounding errors
end

%generates a matrix that contains the possible Fock states in rows; the number in row n is the number of photons in resonator n
%this function needs to be called before anything else
%for example, for the 2 photon, 4 resonator case, resonator_states gets populated with
%     0   0   0   0
%     1   0   0   0                                                                                                                             
%     0   1   0   0                                                                                                                             
%     0   0   1   0                                                                                                                             
%     0   0   0   1                                                                                                                             
%     2   0   0   0                                                                                                                             
%     1   1   0   0                                                                                                                             
%     1   0   1   0                                                                                                                             
%     1   0   0   1                                                                                                                             
%     0   2   0   0                                                                                                                             
%     0   1   1   0
%     0   1   0   1
%     0   0   2   0
%     0   0   1   1
%     0   0   0   2
function generate_resonator_states
  global N_ph;
  global N_res;
  global N_resonator_states;
  global resonator_states;

  N_resonator_states = 1;
  resonator_states = zeros(1,N_res);

  for ph = 1:N_ph %go over all numbers of photons except for zero photon case, which is already covered
    N_resonator_states += nchoosek(N_res+ph-1,ph);
    for st = 1:size(resonator_states)(1) %go over all the states we've already enumerated
      for res = 1:N_res %go over all the resonators
	%come up with a new state by adding a photon to the 'res'th resonator of a previously listed state
	temp = resonator_states(st,:);
	temp(res) += 1; 
	%however, this state may already have made it in to the list

	already_there = 0; %check to see if this state is already in the list
	for k = 1:size(resonator_states)(1)
	  if isequal(resonator_states(k,:),temp)
	    already_there = 1;
	  end
	end

	if already_there == 0 %add it if it's not
	  resonator_states = [resonator_states; temp];
	end
      end
    end
  end
end

function set_resonator_states(new_N_ph, new_N_res)
  global N_ph; %consider up to this many photons
  global N_res; %number of resonators to consider

  if (N_ph != new_N_ph) || (N_res != new_N_res)
    N_ph = new_N_ph;
    N_res = new_N_res;
    generate_resonator_states;
  end
end

% takes a fock state like [1,0,0,1] (one photon in first and last resonators) and returns the equivilent state: [0;0;0;0;0;0;0;0;1;0;0;0;0;0;0]
% only works for the values in the array of states, not superpositions of them
function s = fock_state_to_resonator_state(fock_state)
  global N_resonator_states;
  global resonator_states;
  s = zeros(N_resonator_states, 1);
  s(resonator_state_index(fock_state)) += 1;
end

%reverse of previous function; will work on superpositions of states
function EV = resonator_state_to_EV(state)
  global N_res;
  for n = 1:N_res
    EV(n) = state_expectation_value(state, resonator_number(n));
  end
end

%returns the row number for 'fock_state' in 'resonator_states' matrix
function ind = resonator_state_index(fock_state)
  global resonator_states;
  for n = 1:size(resonator_states)(1)
    if isequal(resonator_states(n,:),fock_state)
      ind = n;
    end
  end
end


% takes a fock_state like [1,1,0,0] (first two QBs in excited state, second two in ground state) and returns the equivilent state
% doesn't work on superpositions of fock_states not superpositions of them
function s = fock_state_to_QB_state(fock_state)
  global N_res; %that's how many QBs we have

  s = [1];

  for n = 1:N_res
    if fock_state(n) == 1
      s = kron(s,[0;1]); %'n'th resonator in excited state
    else
      s = kron(s,[1;0]); %'n'th resonator in gnd state
    end
  end
end

%reverse of previous function; will work on fractional states
function s = QB_state_to_EV(state)
  global N_res;
  for n = 1:N_res
    EV(n) = state_expectation_value(state, QB_number(n));
  end
end

%gives the state of the whole system from fock_states for resonators and QBs
function s = fock_to_state(res_fock,QB_fock)
  s = kron(fock_state_to_resonator_state(res_fock),fock_state_to_QB_state(QB_fock));
end

%gives the density matrix of the whole system from fock_states for resonators and QBs
function dm = fock_to_DM(res_fock,QB_fock)
  dm = diag(fock_to_state(res_fock,QB_fock));
end


%gives the expected number of photons in each resonator for a state or states
function EVs = state_to_resonator_EV(states)
  global N_res;

  for n = 1:N_res %generate all the number operators we'll need
    N_op{n} = kron(resonator_number(n),eye(2^N_res));
  end

  EVs = [];
  for state = states %iterate over the states
    for n = 1:N_res %find the denisty for the state we're on
      EV(n) = state_expectation_value(state, N_op{n});
    end
    EVs = [EVs; EV];%add the denisty we found to the matrix of all EVs
  end
end

%gives the expected number of photons in each resonator for a state or states
function EVs = state_to_QB_EV(states)
  global N_resonator_states;
  global N_res;

  for n = 1:N_res
    N_op{n} = kron(eye(N_resonator_states),QB_number(n));
  end

  EVs = [];
  for state = states %iterate over the states
    for n = 1:N_res %find the denisty for the state we're on
      EV(n) = state_expectation_value(state, N_op{n});
    end
    EVs = [EVs; EV];%add the denisty we found to the matrix of all EVs
  end
end

%next two functions are helper functions; they return the fock state resulting from adding or removing a photon form a given resonator
function st = inc_resonator_fock_state(state, res_index)
  state(res_index) += 1;
  st = state;
end

function st = dec_resonator_fock_state(state, res_index)
  state(res_index) -= 1;
  st = state;
end

%HAD ANNIHILATION AND CREATION BACKWARDS; FIXED NOW, BUT FIGURE OUT WHY THIS WORKS AND CHANGE COMMENTS TO FIT
function m = resonator_annihilation(res_index)
  global N_resonator_states;
  global N_res;
  global N_ph;
  global resonator_states;

  m = zeros(N_resonator_states, N_resonator_states);

  for state = 1:(N_resonator_states-nchoosek(N_res+N_ph-1,N_ph)) %only consider the states with less than N_ph photons
    m(state,resonator_state_index(inc_resonator_fock_state(resonator_states(state,:),res_index))) = sqrt(1+resonator_states(state,:)(res_index));
  end

end

function m = resonator_creation(res_index)
  global N_resonator_states;
  global N_res;
  global N_ph;
  global resonator_states;

  m = zeros(N_resonator_states, N_resonator_states);

  for state = 2:(N_resonator_states) %only consider the states with one or more photons
    if resonator_states(state,:)(res_index) > 0 %can't lower from a state that has nothing in it
      m(state,resonator_state_index(dec_resonator_fock_state(resonator_states(state,:),res_index))) = sqrt(resonator_states(state,:)(res_index));
    end
  end

end

%number operator on 'res_index'th resonator
function m = resonator_number(res_index)
  m = resonator_creation(res_index)*resonator_annihilation(res_index);
end

%the creation operator for the 'QB_index'th QB
%the idea is that we just kroneker product the 2x2 creation matrix for the 'QB_index'th QB with identity QBs for the others
function m = QB_creation(QB_index)
  global N_res; %that's how many QBs we have

  m = [1];

  for n = 1:N_res
    if n == QB_index
      m = kron(m,[0,0;1,0]);
    else
      m = kron(m,eye(2));
    end
  end
end

function m = QB_annihilation(QB_index)
  global N_res; %that's how many QBs we have

  m = [1];

  for n = 1:N_res
    if n == QB_index
      m = kron(m,[0,1;0,0]);
    else
      m = kron(m,eye(2));
    end
  end
end

%number operator on 'QB_index'th QB
function m = QB_number(QB_index)
  m = QB_creation(QB_index)*QB_annihilation(QB_index);
end

%sum of individual QB_number operators
function m = total_QB_number
  global N_res;
  m = 0;

  for n = 1:N_res
    m += QB_creation(n)*QB_annihilation(n);
  end
end

%hamiltonian for just the resonator part of the system
function H = resonator_hamiltonian
  global N_ph;
  global N_res;
  global N_resonator_states;
  global H;
  global h_bar;
  global w0;
  global Wcoup;

  H = zeros(N_resonator_states, N_resonator_states);

  for n = 1:N_res
    H += h_bar*w0*resonator_creation(n)*resonator_annihilation(n); %yeilds the diagonal elements
    %next two 'H +=' lines add in the coupling between the resonators -- the off diagonal elements
    if n > 1
      H += h_bar*Wcoup*resonator_creation(n-1)*resonator_annihilation(n);
    end
    if n < N_res
      H += h_bar*Wcoup*resonator_creation(n+1)*resonator_annihilation(n);
    end
  end
end

%hamiltonian for just the QB part of the system
function H = QB_hamiltonian
  global h_bar;
  global wQB;

  H = h_bar*wQB*total_QB_number;
end

%the hamiltonian for the whole system, with 'interacting' being the lists of QBs coupled to their resonators
%For example, interacting=[1,3,4] would mean QBs 1,3, and 4 are coupled to their resonators
function H = hamiltonian(interacting=[])
  global N_ph;
  global N_res;
  global h_bar;
  global w0;
  global wQB;
  global wInt;
  global Wcoup;
  global N_resonator_states;
  global resonator_states;

  H_res = resonator_hamiltonian;
  H_QB  = QB_hamiltonian;
  I_res = eye(N_resonator_states);
  I_qb  = eye(2^N_res);

  H = kron(H_res, I_qb) + kron(I_res,H_QB);%this is the hamiltonain not including interaction between QBs and resonators

  %now, add in the requested interaction(s)
  for n = 1:size(interacting)(2)
    index = interacting(n);
    H += h_bar*wInt*(kron(resonator_annihilation(index), QB_creation(index)) + kron(resonator_creation(index), QB_annihilation(index)));
  end
end

%evolves 'initial_state' for time 't' with indicies of interacting QBs given in 'interacting'
function psi = evolve_state(initial_state,t,interacting=[])
  global h_bar;

  H = hamiltonian(interacting);

  %method 1, decompose initial_state in to the eigenvectors of 'H', then evolve those -- this seems to be fastest, especially for more than one time value
  [eig_vectors,D] = eig(H);
  eig_values = diag(D);
  psi = zeros(size(H)(1),size(t)(2)); %number of rows equal to number of states, number of columns equal to number of times
  for n = 1:size(eig_vectors)(1)
    projection = eig_vectors(:,n) * (eig_vectors(:,n)'*initial_state); %projection of 'initial_state' on to 'n'th eigenvector
    for m = 1:size(t)(2)
      psi(:,m) += exp(-i*eig_values(n)*t(m)/h_bar)*projection; %evolve weighted 'n'th eigenvector and add to results
    end
  end

  %method 2, exponentiate the matrix -- this seems to be slower; would need to rewrite to handle 't' being an array
%    psi = expm(-i*(H*t)/h_bar)*state;
end;

%evolves 'DM' for time 't' with indicies of interacting QBs given in 'interacting'
function DMs = evolve_DM(DM,times,interacting=[])
  global h_bar;

  H = hamiltonian(interacting);
  DMs = {};

  for t = times
    ev = expm(-i*H*t/h_bar);
    DMs = [DMs, {ev*DM*(ev')}];
  end
end;

%evolves 'initial_state' to times [0,stop_time/(npoints-1),2*stop_time/(npoints-1),...,stop_time] with indicies of
%interacting QBs given in 'interacting' this function works by exponentiating a matrix; that's a high fixed time cost,
%but it only has to be done once because the times are easilly space, so this can be faster than the previous
%function if you're dealing with many evenly spaced time steps
function psi = evolve_state_linspaced_t(initial_state,stop_time,npoints,interacting=[])
  global H0;
  global Hi1;
  global Hi2;
  global h_bar;

  H = hamiltonian(interacting);
  dt = stop_time/(npoints-1);

  ev = expm(-i*H*dt/h_bar);

  psi(:,1) = initial_state;

  for n = 2:npoints
    psi(:,n) = ev*psi(:,n-1);
  end
end;

function DMs = evolve_DM_linspaced_t(initial_DM,stop_time,npoints,interacting=[])
  global h_bar;

  H = hamiltonian(interacting);
  dt = stop_time/(npoints-1);

  evL = expm(-i*H*dt/h_bar);
  evR = evL';

  DMs{1} = initial_DM;

  for n = 2:npoints
    DMs{n} = evL*DMs{n-1}*evR;
  end
end;

%ONLY MULTI-TIME EVOLUTION FUNCTION THAT SUPPORTS RELAXATION!!
function DMs = relax_evolve_DM_linspaced_t(initial_DM,stop_time,npoints,interacting=[],store_intermediate_points=true)
  global h_bar;
  global N_res;
  global N_resonator_states;
  global N_ph;
  global QBRelaxationTime;
  global resRelaxationTime;
  global resonator_states;

  H = hamiltonian(interacting);
  dt = stop_time/(npoints-1);

  evL = expm(-i*H*dt/h_bar);
  evR = evL';
  
  I_res = eye(N_resonator_states);
  I_qb  = eye(2^N_res);

  %store the E0 and E1 operators
  E0s = {};
  E1s = {};

  %make the E0 and E1 Kruas operators for amplitude damping of the qubits (see Nielsen & Chuang eqn 8.108, p 380)
  P_QB_relax = 1-exp(-dt/QBRelaxationTime); 
  for n = 1:N_res
    E0s{n} = [1];
    E1s{n} = [1];
    for m = 1:N_res
      if m == n
	E0s{n} = kron(E0s{n},[1,0;0,sqrt(1-P_QB_relax)]);
	E1s{n} = kron(E1s{n},[0,sqrt(P_QB_relax);0,0]);
      else
	E0s{n} = kron(E0s{n},eye(2));
	E1s{n} = kron(E1s{n},eye(2));
      end
    end
    %turn them in to operators that can act on the full density matrix
    E0s{n} = kron(I_res, E0s{n});
    E1s{n} = kron(I_res, E1s{n});
  end

  %make the E0 and E1 Kruas operators for amplitude damping of the resonators (see Nielsen & Chuang eqn 8.108, p 380)
  for n = 1:N_resonator_states
    for m = 1:N_res
      num_photons = resonator_states(n,m);
      if num_photons > 0
	P_res_relax = 1-exp(-dt*num_photons/resRelaxationTime); %relaxation time decreases like 1/(number of photons)
	%'relaxing_from' is the index of the state we're currently in
	relaxing_from = n;
	%'relaxing_to' is the index of the state we're relaxing to
	relaxing_to = resonator_state_index(dec_resonator_fock_state(resonator_states(relaxing_from,:),m)); %we're relaxing resonator m
	%E0 is all ones on the diagonal except it's sqrt(1-P_res_relax) for the state that's relaxing (that state has index 'n')
	E0_temp = diag(ones(1,N_resonator_states));
	E0_temp(relaxing_from,relaxing_from) = sqrt(1-P_res_relax);
	%E1 is all zeros except it's sqrt(P_res_relax) connecting the state that's relaxing to the state it's relaxing to
	E1_temp = zeros(N_resonator_states, N_resonator_states);
	E1_temp(relaxing_to,relaxing_from) = sqrt(P_res_relax);
	%turn them in to operators that can act on the full density matrix
	E0_temp = kron(E0_temp, I_qb);
	E1_temp = kron(E1_temp, I_qb);
	%put them in the cell array
	E0s = [E0s, {E0_temp}];
	E1s = [E1s, {E1_temp}];
      end
    end
  end

  if store_intermediate_points
    DMs{1} = initial_DM;
    for n = 2:npoints
      DMs{n} = evL*DMs{n-1}*evR; %evolve one timestep
      for m = 1:size(E0s)(2) %loop over all the relaxation operators
	DMs{n} = E0s{m}*DMs{n}*(E0s{m}') + E1s{m}*DMs{n}*(E1s{m}');
      end
    end
  else
    DMs = initial_DM;
    for n = 2:npoints
      DMs = evL*DMs*evR; %evolve one timestep
      for m = 1:size(E0s)(2) %loop over all the relaxation operators
	DMs = E0s{m}*DMs*(E0s{m}') + E1s{m}*DMs*(E1s{m}');
      end
    end
  end
end;

%only evolves for one time
function DM = relax_evolve_DM(initial_DM,stop_time,interacting=[])
  global RelEvlStepSz;
  
  if stop_time > 0
    npoints = ceil(stop_time/RelEvlStepSz);
    DM = relax_evolve_DM_linspaced_t(initial_DM,stop_time,npoints,interacting,false);
  else %handles case where stop_time == 0
    DM = initial_DM;
  end
end

%takes in a state or array of states and returns an array of probabilities that those states are in the given fock state
function ps = prob_state_in_fock_state(states, resonator_fock_state, QB_fock_state)
  ps = [];
  fock_state = EV_to_state(resonator_fock_state,QB_fock_state);
  for state = states
    ps = [ps , abs(state'*fock_state)^2];
  end
end

%takes in a density matrix or cell array of density matrices and returns an array of probabilities that those states are in the given fock state  
function ps = prob_DM_in_fock_state(DMs, resonator_fock_state, QB_fock_state)
  ps = [];
  fock_state_operator = fock_to_DM(resonator_fock_state, QB_fock_state);
  for DM = DMs
    ps = [ps , DM_expectation_value(DM{1}, fock_state_operator)];
  end
end

%%%%%%%%%%%%%%%%%%% SIMULATIONS OF EXPERIMENTS %%%%%%%%%%%%%%%%%%%


%starts with both QBs excited, swaps them in to the resonators, waits the given time, and then swaps out for swap time/sqrt(2)
%this simulates the Hong–Ou–Mandel experiment for two resonators and two qubits WITHOUT RELAXATION
%it returns g2, the joint switching probabilities divided by the individual switching probabilities (P11/P1A/P1B)
%'wait_times' is the array of times you want g2 at
function g2 = HOM(wait_times)
  global h_bar;
  global wInt;

  set_resonator_states(2, 2);

  QB_swap_t = 2*pi*1/wInt/4;

  initial_state = EV_to_state([0,0],[1,1]);
  swap_in = expm(-i*hamiltonian([1,2])*QB_swap_t/h_bar);
  swap_out = expm(-i*hamiltonian([1,2])*QB_swap_t/sqrt(2)/h_bar); %two photons in the resonators mean that swap time is multiplied by 1/sqrt(2)

  final_states = [];
  for wait_time = wait_times
    wait = expm(-i*hamiltonian()*wait_time/h_bar);
    final_states = [final_states, swap_out*wait*swap_in * initial_state];
  end

  %probability junction A is excited
  P1A = prob_state_in_fock_state(final_states, [0,0], [1,1]) + prob_state_in_fock_state(final_states, [1,0], [0,1]) + prob_state_in_fock_state(final_states, [0,1], [0,1]);
  %probability junction B is excited
  P1B = prob_state_in_fock_state(final_states, [0,0], [1,1]) + prob_state_in_fock_state(final_states, [1,0], [1,0]) + prob_state_in_fock_state(final_states, [0,1], [1,0]);
  %probability both junctions are excited
  P11 = prob_state_in_fock_state(final_states, [0,0], [1,1]);
  g2 = P11./P1A./P1B;
end

%starts with both QBs excited, swaps them in to the resonators, waits the given time, and then swaps out for swap time/sqrt(2)
%this simulates the Hong–Ou–Mandel experiment for two resonators and two qubits WITH RELAXATION
%it returns g2, the joint switching probabilities divided by the individual switching probabilities (P11/P1A/P1B)
%'wait_times' is the array of times you want g2 at
function g2 = HOMrelax(wait_length)
  global wInt;

  set_resonator_states(2, 2);

  QB_swap_t = 2*pi*1/wInt/4;

  initial_DM = fock_to_DM([0,0],[1,1]);
  after_swap_in = relax_evolve_DM_linspaced_t(initial_DM,QB_swap_t,QB_swap_t/10e-12,[1,2]){end};
  after_wait_times = relax_evolve_DM_linspaced_t(after_swap_in,wait_length,wait_length/10e-12);
  after_wait_times = after_wait_times(1:500:end); %throw out 499/500 points TODO: DO THIS MORE CLEANLY; ADD OPTION

  after_swapping_out = {};
  for aw = after_wait_times
    after_swapping_out = [after_swapping_out, {relax_evolve_DM_linspaced_t(aw{1},QB_swap_t/sqrt(2),QB_swap_t/sqrt(2)/10e-12,[1,2]){end}}];
  end


  %probability junction A is excited
  P1A = prob_DM_in_fock_state(after_swapping_out, [0,0], [1,1]) + prob_DM_in_fock_state(after_swapping_out, [1,0], [0,1]) + prob_DM_in_fock_state(after_swapping_out, [0,1], [0,1]);
  %probability junction B is excited
  P1B = prob_DM_in_fock_state(after_swapping_out, [0,0], [1,1]) + prob_DM_in_fock_state(after_swapping_out, [1,0], [1,0]) + prob_DM_in_fock_state(after_swapping_out, [0,1], [1,0]);
  %probability both junctions are excited
  P11 = prob_DM_in_fock_state(after_swapping_out, [0,0], [1,1]);
  g2 = P11./P1A./P1B;
end

%TODO: WRITE EXPLANATION
function g2 = HOMrelax_delay
  global wInt;
  global RelEvlStepSz;

  set_resonator_states(2, 2);

  QB1_start_t_min = 10e-9;
  QB1_start_t_max = 70e-9;
  QB2_start_t     = 40e-9;
  max_delay       = 10e-9;

  QB_swap_in_t  = 2*pi*1/wInt/4;
  QB_swap_out_t = 2*pi*1/wInt/4/sqrt(2);

  QB1_start_times = linspace(QB1_start_t_min,QB1_start_t_max,101); %might run in to problems if we start right at zero
  QB2_rel_start_t = QB2_start_t - QB1_start_t_min; %relative times are times after QB1_start_t_min
  QB2_rel_end_t   = QB2_rel_start_t + QB_swap_out_t;

  initial_DM = fock_to_DM([0,0],[1,1]);
  after_swap_in  = relax_evolve_DM(initial_DM,    QB_swap_in_t,    [1,2]);
  after_min_time = relax_evolve_DM(after_swap_in, QB1_start_t_min, []);

  after_swapping_out = {};
n = 1
  for QB1_start_t = QB1_start_times
n = n + 1
    QB1_rel_start_t = QB1_start_t - QB1_start_t_min;
    QB1_rel_end_t   = QB1_rel_start_t + QB_swap_out_t;

    %by symetry, it doesn't really matter which one starts first or second, only that one starts first or second
    %so QBa will be the one starting first, and QBb will be the one starting second
    QBa_s = min(QB1_rel_start_t, QB2_rel_start_t);
    QBa_e = min(QB1_rel_end_t, QB2_rel_end_t);
    QBb_s = max(QB1_rel_start_t, QB2_rel_start_t);
    QBb_e = max(QB1_rel_end_t, QB2_rel_end_t);

    temp = relax_evolve_DM(after_min_time, QBa_s, []); %evolve up to when QBa starts interacting
    if QBb_s < QBa_e %QBb starts before QBa ends
      temp = relax_evolve_DM(temp, QBb_s - QBa_s, [1]); %evolve up to when QBb starts interacting
      temp = relax_evolve_DM(temp, QBa_e - QBb_s, [1,2]); %evolve up to when QBa stops interacting
      temp = relax_evolve_DM(temp, QBb_e - QBa_e, [2]); %evolve up to when QBb stops interacting
    else %QBb starts after QBa ends
      temp = relax_evolve_DM(temp, QBa_e - QBa_s, [1]); %evolve up to when QBa stops interacting
      temp = relax_evolve_DM(temp, QBb_s - QBa_e, []); %evolve up to when QBb starts interacting
      temp = relax_evolve_DM(temp, QBb_e - QBb_s, [2]); %evolve up to when QBb stops interacting
    end
%      QB1_start_t = delay_t + QB1_start_t_min;
%      if delay_t + QB1_start_t_min + QB_swap_out_t < QB2_start_t %QBs never interact at the same time; QB1 is done before QB2 starts
%        temp = relax_evolve_DM(after_min_time, delay_t, []);
%        temp = relax_evolve_DM(temp,           QB_swap_out_t, [1]);
%        temp = relax_evolve_DM(temp,           QB2_start_t - QB_swap_out_t - QB1_start_t, []);
%        temp = relax_evolve_DM(temp,           QB_swap_out_t, [2]);
%      elseif delay_t + QB1_start_t_min > QB2_start_t + QB_swap_out_t %QBs never interact at the same time; QB2 is done before QB1 starts
%        temp = relax_evolve_DM(after_min_time, QB2_start_t - QB1_start_t_min, []);
%        temp = relax_evolve_DM(temp,           QB_swap_out_t, [2]);
%        temp = relax_evolve_DM(temp,           QB1_start_t - QB2_start_t - QB_swap_out_t, []);
%        temp = relax_evolve_DM(temp,           QB_swap_out_t, [1]);
%      elseif delay_t + QB1_start_t_min < QB2_start_t %QBs' interaction periods overlap; QB1 starts first
%        temp = relax_evolve_DM(after_min_time, delay_t, []);
%        temp = relax_evolve_DM(temp,           QB2_start_t - QB1_start_t, [1]);
%        temp = relax_evolve_DM(temp,           QB1_start_t + QB_swap_out_t - QB2_start_t, [1,2]);
%        temp = relax_evolve_DM(temp,           QB1_start_t + QB_swap_out_t - QB2_start_t, [1,2]);
%      else %QBs' interaction periods overlap; QB2 starts first
%        temp = relax_evolve_DM(after_min_time, QB2_start_t - QB1_start_t_min, []);
%        temp = relax_evolve_DM(temp,           QB1_start_t - QB2_start_t, [2]);
%        temp = relax_evolve_DM(temp,           QB1_start_t - QB2_start_t, [2]);
%  NOT FINISHED
%      end

    after_swapping_out = [after_swapping_out, {temp}];
  end


  %probability junction A is excited
  P1A = prob_DM_in_fock_state(after_swapping_out, [0,0], [1,1]) + prob_DM_in_fock_state(after_swapping_out, [1,0], [0,1]) + prob_DM_in_fock_state(after_swapping_out, [0,1], [0,1]);
  %probability junction B is excited
  P1B = prob_DM_in_fock_state(after_swapping_out, [0,0], [1,1]) + prob_DM_in_fock_state(after_swapping_out, [1,0], [1,0]) + prob_DM_in_fock_state(after_swapping_out, [0,1], [1,0]);
  %probability both junctions are excited
  P11 = prob_DM_in_fock_state(after_swapping_out, [0,0], [1,1]);
  g2 = P11./P1A./P1B;
end

%%%%%%%%%%%%%%%%%%% TEMPORARY FUNCTIONS %%%%%%%%%%%%%%%%%%%


function p11 = p11s(wait_length)
  global wInt;

  set_resonator_states(2, 2);

  QB_swap_t = 2*pi*1/wInt/4;

  initial_DM = fock_to_DM([0,0],[1,1]);
  after_swap_in = relax_evolve_DM_linspaced_t(initial_DM,QB_swap_t,QB_swap_t/10e-12,[1,2]){end};
  after_wait_times = relax_evolve_DM_linspaced_t(after_swap_in,wait_length,wait_length/10e-12);
  after_wait_times = after_wait_times(3500:50:end); %throw out 499/500 points TODO: DO THIS MORE CLEANLY; ADD OPTION

  after_swapping_out = {};
  for aw = after_wait_times
    after_swapping_out = [after_swapping_out, {relax_evolve_DM_linspaced_t(aw{1},QB_swap_t/sqrt(2),QB_swap_t/sqrt(2)/10e-12,[1,2]){end}}];
  end

  p11 = prob_DM_in_fock_state(after_swapping_out, [0,0], [1,1]);
end