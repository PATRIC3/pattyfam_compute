TOP_DIR = ../..
include $(TOP_DIR)/tools/Makefile.common

DEPLOY_RUNTIME ?= /kb/runtime
TARGET ?= /kb/deployment

APP_CXX = propagate_names
BIN_CXX = $(addprefix $(BIN_DIR)/,$(APP_CXX))
DEPLOY_CXX = $(addprefix $(TARGET)/bin,$(APP_CXX))

SRC_PERL = $(wildcard scripts/*.pl)
BIN_PERL = $(addprefix $(BIN_DIR)/,$(basename $(notdir $(SRC_PERL))))
DEPLOY_PERL = $(addprefix $(TARGET)/bin/,$(basename $(notdir $(SRC_PERL))))

SRC_SERVICE_PERL = $(wildcard service-scripts/*.pl)
BIN_SERVICE_PERL = $(addprefix $(BIN_DIR)/,$(basename $(notdir $(SRC_SERVICE_PERL))))
DEPLOY_SERVICE_PERL = $(addprefix $(SERVICE_DIR)/bin/,$(basename $(notdir $(SRC_SERVICE_PERL))))

CLIENT_TESTS = $(wildcard t/client-tests/*.t)
SERVER_TESTS = $(wildcard t/server-tests/*.t)
PROD_TESTS = $(wildcard t/prod-tests/*.t)

STARMAN_WORKERS = 8
STARMAN_MAX_REQUESTS = 100

TPAGE_ARGS = --define kb_top=$(TARGET) --define kb_runtime=$(DEPLOY_RUNTIME) --define kb_service_name=$(SERVICE) \
	--define kb_service_port=$(SERVICE_PORT) --define kb_service_dir=$(SERVICE_DIR) \
	--define kb_sphinx_port=$(SPHINX_PORT) --define kb_sphinx_host=$(SPHINX_HOST) \
	--define kb_starman_workers=$(STARMAN_WORKERS) \
	--define kb_starman_max_requests=$(STARMAN_MAX_REQUESTS)

all: bin 

bin: $(BIN_PERL) $(BIN_SERVICE_PERL) $(BIN_CXX)

#PROFILE = -pg
OPT = -O3
DEBUG = -g
INC = $(BOOST_INC) $(TBB_FLAGS) $(NUDB_INCLUDE) $(CMPH_INCLUDE) -I../signature_kmers/src
CXXFLAGS = $(PROFILE) $(DEBUG) $(OPT) $(INC)
LDFLAGS = -Wl,-rpath,$(BOOST)/lib -Wl,-rpath,$(CMPH)/lib $(PROFILE)

LIBS = $(BOOST_LIBS) $(TBB_LIBS) $(CMPH_LIB) $(LOG4CXX_LIB)

LOG4CXX_LIB = -llog4cxx

BOOST = $(KB_RUNTIME)/boost-latest

BOOST_INC = -I$(BOOST)/include
BOOST_LIBS = \
	-L $(BOOST)/lib \
	-lboost_program_options \
	-lboost_filesystem \
	-lboost_thread \
	-lboost_regex
#
# Parallel algorithms on ubuntu result in TBB deprecation messages
#
TBB_FLAGS = -DTBB_SUPPRESS_DEPRECATED_MESSAGES=1 -DTBB_PREVIEW_CONCURRENT_ORDERED_CONTAINERS=1
TBB_LIBS = -ltbbmalloc -ltbb

PROPAGATE_NAMES_OBJS = src/propagate_names.o
propagate_names: $(PROPAGATE_NAMES_OBJS)
	$(CXX) $(LDFLAGS) -o $@ $(PROPAGATE_NAMES_OBJS) $(LIBS)


deploy: deploy-all
deploy-all: deploy-client 
deploy-client: deploy-libs deploy-scripts deploy-docs

deploy-dir:
	if [ ! -d $(SERVICE_DIR) ] ; then mkdir $(SERVICE_DIR) ; fi
	if [ ! -d $(SERVICE_DIR)/webroot ] ; then mkdir $(SERVICE_DIR)/webroot ; fi
	if [ ! -d $(SERVICE_DIR)/bin ] ; then mkdir $(SERVICE_DIR)/bin ; fi
	if [ ! -d $(SERVICE_DIR)/sphinx ] ; then mkdir $(SERVICE_DIR)/sphinx ; fi

deploy-docs: 

test: test-client
	echo "running client and script tests"

# What does it mean to test a client. This is a test of a client
# library. If it is a client-server module, then it should be
# run against a running server. You can say that this also tests
# the server, and I agree. You can add a test-server dependancy
# to the test-client target if it makes sense to you. This test
# example assumes there is already a tested running server.
test-client:
	# run each client test
	for t in $(CLIENT_TESTS) ; do \
		if [ -f $$t ] ; then \
			$(DEPLOY_RUNTIME)/bin/perl $$t ; \
			if [ $$? -ne 0 ] ; then \
				exit 1 ; \
			fi \
		fi \
	done

test-server:
	# run each server test
	for t in $(SERVER_TESTS) ; do \
		if [ -f $$t ] ; then \
			echo Running $$t ; \
			$(DEPLOY_RUNTIME)/bin/perl $$t ; \
			if [ $$? -ne 0 ] ; then \
				exit 1 ; \
			fi \
		fi \
	done

test-prod-server:
	# run each prod test
	for t in $(PROD_TESTS) ; do \
		if [ -f $$t ] ; then \
			echo Running $$t ; \
			$(DEPLOY_RUNTIME)/bin/perl $$t ; \
			if [ $$? -ne 0 ] ; then \
				exit 1 ; \
			fi \
		fi \
	done

clean:
	ant clean

include $(TOP_DIR)/tools/Makefile.common.rules
