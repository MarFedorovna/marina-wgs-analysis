import datetime as dt
from pathlib import Path
import signal
import sys
import time
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import continue_packet_views_fast_v2_1 as controller


@unittest.skipUnless(sys.platform=='linux','POSIX whole-parent deadline and address-space guards')
class ParentLimits(unittest.TestCase):
    def test_parent_merge_wait_is_interrupted_and_timer_restored(self):
        previous=signal.getsignal(signal.SIGALRM)
        start=time.monotonic()
        with self.assertRaisesRegex(TimeoutError,'parent deadline'):
            with controller.parent_limits(controller.now()+dt.timedelta(seconds=.05)):
                time.sleep(2)
        self.assertLess(time.monotonic()-start,1)
        self.assertEqual(signal.getsignal(signal.SIGALRM),previous)
        self.assertEqual(signal.getitimer(signal.ITIMER_REAL),(0.,0.))

    def test_two_gib_parent_limit_restored_after_success_and_exception(self):
        import resource
        previous=resource.getrlimit(resource.RLIMIT_AS)
        with controller.parent_limits(controller.now()+dt.timedelta(seconds=2)):
            self.assertLessEqual(resource.getrlimit(resource.RLIMIT_AS)[0],2*1024**3)
        self.assertEqual(resource.getrlimit(resource.RLIMIT_AS),previous)
        with self.assertRaises(RuntimeError):
            with controller.parent_limits(controller.now()+dt.timedelta(seconds=2)):raise RuntimeError('test only')
        self.assertEqual(resource.getrlimit(resource.RLIMIT_AS),previous)
        self.assertEqual(signal.getitimer(signal.ITIMER_REAL),(0.,0.))

    def test_expired_parent_never_runs_body(self):
        with self.assertRaisesRegex(ValueError,'parent deadline'):
            with controller.parent_limits(controller.now()-dt.timedelta(seconds=1)):self.fail('Must never run')


if __name__=='__main__':unittest.main()
