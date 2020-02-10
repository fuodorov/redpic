void loadbeam( struct particle*, struct parametrs);

int run_12(struct particle*, struct parametrs, struct object*);

int run_0 (struct particle*, struct parametrs, struct object*);

struct gauss_bunch loadbunch(void);

double gauss(double);

void save_data(void);

struct parametrs par_to_in(struct parametrs, struct object* );

struct object*    obj_to_in(struct parametrs, struct object* );

int where_is_particle(struct particle, struct parametrs, struct object*);

struct vector E_obj(int, struct particle, struct parametrs, struct object);

struct vector H_obj(int, struct particle, struct parametrs, struct object);

void draw_beam(struct object*);
