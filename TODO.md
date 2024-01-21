# Remaining TODO

General:
* notes scattered in comments and docstrings (sorry!)
* Docker-ification. More or less needs to start from scratch (there should be existing images for SeqRepo and UTA, not sure how up-to-date they are).
* Add extra stuff that appears in mapping JSON objects.
* Currently using VRS 2.0a-based libraries. For lifting back to VRS 1.3, some basic post-processing should be fine (annoying but shouldn't be too trivial)

Alignment:
* Pretty sure this is mostly done. Haven't tested exhaustively, though.
* Need to sufficiently mock/patch things in tests

Transcript selection:
* IndexError in calculating offset on lots of new (2023) scoresets.
* Tests will need some extensive mocking (or cassettes?) for reliance on UTA and other external dependencies

VRS mapping:
* In general, this stuff is still pretty rough
* Finish the SeqRepo storage workaround
* A fair amount of small questions about conditions written to handle specific scoresets/edge cases
* More testing. Can be ready for CI by mocking the SequenceStore class (or using SeqRepoRESTDataProxy and cassettes).
