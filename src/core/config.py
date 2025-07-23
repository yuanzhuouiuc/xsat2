from dataclasses import dataclass

@dataclass
class SolverConfig:
    niter: int = 100
    method: str = 'powell'
    multi: bool = True
    step_size: float = 10.0