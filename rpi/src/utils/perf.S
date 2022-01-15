.globl get_cycle_count
get_cycle_count:
    mrc p15, 0, r0, c9, c13, 0
    mov pc, lr


